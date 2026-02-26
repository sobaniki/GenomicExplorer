from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, List


def build_runner_command(
    *,
    env: Dict[str, Any],
    language: str,          # "r"
    runner_path: Path,
    params_json_path: Path,
    out_dir: Path,
) -> List[str]:
    mode = (env.get("mode") or "system").lower()
    language = language.lower()

    if language in {"r", "rscript"}:
        base = [str(env.get("rscript") or "Rscript"), str(runner_path)]
    else:
        base = [str(env.get("python") or "python"), str(runner_path)]

    args = ["--params", str(params_json_path), "--out", str(out_dir)]

    if mode == "system":
        return base + args
    if mode == "conda":
        name = env.get("name")
        if not name:
            raise ValueError("env.mode=conda requires env.name")
        return ["conda", "run", "-n", str(name)] + base + args

    # venv 等は必要になったら増やす
    return base + args

