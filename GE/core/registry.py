from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional


@dataclass(frozen=True)
class PluginRecord:
    plugin_id: str
    plugin_dir: Path
    env: Dict[str, Any]  # {"mode": "...", "name": "...", ...}

@dataclass(frozen=True)
class Registry:
    path: Path
    plugins: Dict[str, PluginRecord]

    @staticmethod
    def load(path: str | Path) -> "Registry":
        path = Path(path)
        raw = json.loads(path.read_text(encoding="utf-8"))
        items = raw.get("plugins", [])
        plugins: Dict[str, PluginRecord] = {}
        for it in items:
            pid = it["id"]
            ppath = Path(it["path"]).expanduser()
            pdir = (path.parent / ppath).resolve() if not ppath.is_absolute() else ppath.resolve()
            env = dict(it.get("env", {}) or {})
            plugins[pid] = PluginRecord(plugin_id=pid, plugin_dir=pdir, env=env)
        return Registry(path=path, plugins=plugins)

    def get(self, plugin_id: str) -> PluginRecord:
        if plugin_id not in self.plugins:
            raise KeyError(f"Plugin not found in registry: {plugin_id}")
        return self.plugins[plugin_id]
