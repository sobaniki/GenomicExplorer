from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional
import yaml


@dataclass(frozen=True)
class PluginManifest:
    plugin_dir: Path
    raw: Dict[str, Any]

    @property
    def plugin_id(self) -> str:
        return str(self.raw.get("id", self.plugin_dir.name))

    @property
    def parameters_schema(self) -> Dict[str, Any]:
        return dict(self.raw.get("parameters_schema", {}) or {})

    @property
    def outputs_schema(self) -> Dict[str, Any]:
        return dict(self.raw.get("outputs_schema", {}) or {})

    @property
    def runner(self) -> Dict[str, Any]:
        return dict(self.raw.get("runner", {}) or {})


def load_manifest(plugin_dir: str | Path) -> PluginManifest:
    plugin_dir = Path(plugin_dir)
    yml = plugin_dir / "plugin.yaml"
    if not yml.exists():
        raise FileNotFoundError(f"plugin.yaml not found: {yml}")
    try:
        raw = yaml.safe_load(yml.read_text(encoding="utf-8")) or {}
    except yaml.YAMLError as e:
        # Keep this message strictly single-line string literals to avoid
        # accidental syntax errors from embedded newlines.
        msg = (
            f"Failed to parse plugin.yaml: {yml}\n"
            f"YAML error: {e}\n"
            "Tip: If the name contains a colon (:) or other punctuation, wrap it in double quotes, e.g.\n"
            "  name: \"My plugin: description\"\n"
        )
        raise ValueError(msg) from e
    if not isinstance(raw, dict):
        raise ValueError(f"plugin.yaml must be a mapping: {yml}")
    return PluginManifest(plugin_dir=plugin_dir, raw=raw)
