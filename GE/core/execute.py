from __future__ import annotations

import json
import subprocess
import shutil
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

from .registry import Registry
from .manifest import load_manifest
from .runner_cmd import build_runner_command


class PluginRunError(RuntimeError):
    def __init__(self, message: str, out_dir, cmd, returncode=None):
        super().__init__(message)
        self.out_dir = out_dir
        self.cmd = cmd
        self.returncode = returncode


def _write_text(path: Path, text: str) -> None:
    path.write_text(text or "", encoding="utf-8", errors="ignore")


def run_plugin(
    registry_path: str | Path,
    plugin_id: str,
    parameters: Dict[str, Any],
    work_dir: str | Path,
) -> Path:
    reg = Registry.load(registry_path)
    rec = reg.get(plugin_id)

    manifest = load_manifest(rec.plugin_dir)
    runner = manifest.runner
    language = (runner.get("language") or "r").lower()
    entry = runner.get("entry")
    if not entry:
        raise RuntimeError("plugin.yaml runner.entry is required")

    work_dir = Path(work_dir).resolve()
    out_dir = work_dir / "out"
    out_dir.mkdir(parents=True, exist_ok=True)

    params_json = work_dir / "params.json"
    _write_text(params_json, json.dumps(parameters, ensure_ascii=False, indent=2))

    # Always keep params alongside outputs for downstream integration (Compare/Integrate).
    _write_text(out_dir / "params.json", json.dumps(parameters, ensure_ascii=False, indent=2))

    runner_path = (rec.plugin_dir / entry).resolve()
    if not runner_path.exists():
        raise RuntimeError(f"Runner not found: {runner_path}")

    cmd = build_runner_command(
        env=rec.env,
        language=language,
        runner_path=runner_path,
        params_json_path=params_json,
        out_dir=out_dir,
    )

    # 事前に実行情報を保存
    _write_text(out_dir / "cmd.txt", " ".join(cmd))
    _write_text(out_dir / "plugin_dir.txt", str(rec.plugin_dir))
    _write_text(out_dir / "work_dir.txt", str(work_dir))

    # Run metadata (for downstream integration / provenance)
    try:
        from datetime import datetime, timezone
        run_meta = {
            "plugin_id": str(plugin_id),
            "plugin_name": str(manifest.raw.get("name", "")),
            "plugin_version": str(manifest.raw.get("version", "")),
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "env": str(rec.env),
            "runner_language": str(language),
            "runner_entry": str(entry),
            "registry_path": str(Path(registry_path).resolve()),
            "plugin_dir": str(rec.plugin_dir),
            "work_dir": str(work_dir),
        }
        _write_text(out_dir / "run_meta.json", json.dumps(run_meta, ensure_ascii=False, indent=2))
    except Exception as _e:
        # Never fail the run just because metadata could not be written
        _write_text(out_dir / "run_meta_error.txt", str(_e))

    try:
        # GUI から呼ばれるケースでは、subprocess.run() をメインスレッドで待つと
        # Qt のイベントループが止まり「応答なし…」になりやすい。
        # ここでは *同期APIのまま* Qt(QProcess + QEventLoop) で待機できる場合は
        # それを優先し、CUI 実行など Qt が無い環境では従来どおり subprocess.run を使う。

        def _run_with_qprocess() -> tuple[int, str, str]:
            try:
                # PySide6 (GUI) / PyQt5 (optional) いずれでも動くように遅延 import
                try:
                    from PySide6.QtCore import QEventLoop, QProcess  # type: ignore
                except Exception:
                    from PyQt5.QtCore import QEventLoop, QProcess  # type: ignore
            except Exception:
                raise ImportError

            p = QProcess()
            p.setWorkingDirectory(str(rec.plugin_dir))
            p.setProcessChannelMode(QProcess.SeparateChannels)

            program = cmd[0]
            args = cmd[1:]

            loop = QEventLoop()
            started_ok = {"ok": False}
            start_error = {"msg": ""}

            def _on_started():
                started_ok["ok"] = True

            def _on_error(_):
                # start失敗など
                try:
                    start_error["msg"] = str(p.errorString())
                except Exception:
                    start_error["msg"] = "Failed to start process"
                loop.quit()

            def _on_finished(_code, _status):
                loop.quit()

            p.started.connect(_on_started)
            p.errorOccurred.connect(_on_error)
            p.finished.connect(_on_finished)

            p.start(program, args)
            loop.exec_()  # GUIイベントを回しつつ待つ

            if not started_ok["ok"] and start_error["msg"]:
                raise RuntimeError(start_error["msg"])

            # finished 後にまとめて取得
            stdout = bytes(p.readAllStandardOutput()).decode("utf-8", errors="replace")
            stderr = bytes(p.readAllStandardError()).decode("utf-8", errors="replace")
            rc = int(p.exitCode())
            return rc, stdout, stderr

        try:
            rc, out_s, err_s = _run_with_qprocess()
            class _Proc:
                returncode = rc
                stdout = out_s
                stderr = err_s
            proc = _Proc()
        except ImportError:
            # Qt が無い(=CUI) / importできない場合は従来どおり
            proc = subprocess.run(
                cmd,
                cwd=str(rec.plugin_dir),
                capture_output=True,
                text=True,
                check=False,
            )
    except Exception as e:
        # conda が PATH にない等、起動前に落ちるケース
        _write_text(out_dir / "stdout.txt", "")
        _write_text(out_dir / "stderr.txt", str(e))
        _write_text(out_dir / "error.txt", f"Failed to start process: {e}\ncmd={' '.join(cmd)}\n")
        raise PluginRunError(
            f"Failed to start runner: {e}",
            out_dir=out_dir,
            cmd=cmd,
            returncode=None,
        ) from e

    _write_text(out_dir / "stdout.txt", proc.stdout or "")
    _write_text(out_dir / "stderr.txt", proc.stderr or "")

    if proc.returncode != 0:
        # 失敗時も見やすい要約を例外に含める
        tail = (proc.stderr or proc.stdout or "")[-1500:]
        _write_text(out_dir / "error.txt",
                    f"Runner failed (code={proc.returncode})\ncmd={' '.join(cmd)}\n\n--- tail ---\n{tail}\n")
        raise PluginRunError(
            f"Runner failed (code={proc.returncode}). See out/stderr.txt (and error.txt). Tail:\n{tail}",
            out_dir=out_dir,
            cmd=cmd,
            returncode=proc.returncode,
        )

    # --- Optional export ---
    # Many GUI panels expose "Export folder (optional)". Historically this field was
    # passed to runners but not applied consistently. To make behavior uniform across
    # all plugins (R and Python runners), we implement export here at the core layer.
    #
    # Supported parameter keys (best-effort): export_dir / export_folder / export_path.
    export_dir_raw = (
        parameters.get("export_dir")
        or parameters.get("export_folder")
        or parameters.get("export_path")
        or ""
    )
    export_dir_raw = str(export_dir_raw).strip()
    if export_dir_raw:
        export_root = Path(export_dir_raw).expanduser().resolve()
        export_root.mkdir(parents=True, exist_ok=True)
        if not export_root.is_dir():
            raise RuntimeError(f"Export folder is not a directory: {export_root}")

        stamp = time.strftime("%Y%m%d_%H%M%S")
        dest = export_root / f"{plugin_id}_{stamp}"
        # If the same second is used, add a suffix.
        if dest.exists():
            i = 2
            while (export_root / f"{plugin_id}_{stamp}_{i}").exists():
                i += 1
            dest = export_root / f"{plugin_id}_{stamp}_{i}"

        shutil.copytree(out_dir, dest)
        # Keep a pointer to the original temp location for debugging.
        _write_text(dest / "_source_work_dir.txt", str(work_dir))
        return dest

    return out_dir

