"""YAML configuration loader.

Loads the simulation configuration plus detector parameters.
Detector-specific fields are placeholders that will be filled in later.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

DEFAULTS_DIR = Path(__file__).resolve().parent.parent.parent.parent / "config"


def _deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge ``override`` into ``base`` (shallow-cloned)."""
    out = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _deep_merge(out[key], value)
        else:
            out[key] = value
    return out


class Config:
    """Simple attribute-style access wrapper over a dict."""

    def __init__(self, data: dict):
        self._data = data

    def get(self, path: str, default=None):
        node = self._data
        for part in path.split("."):
            if not isinstance(node, dict) or part not in node:
                return default
            node = node[part]
        return node

    def as_dict(self) -> dict:
        return self._data


def load_config(path: str | Path | None = None) -> Config:
    """Load a YAML config, merging over the default config (if present)."""
    default_file = DEFAULTS_DIR / "settings.yaml"
    data: dict[str, Any] = {}
    if default_file.exists():
        with open(default_file, "r", encoding="utf-8") as fh:
            data = yaml.safe_load(fh) or {}

    if path is not None:
        cfg_path = Path(path)
        if not cfg_path.exists():
            raise FileNotFoundError(f"Config not found: {cfg_path}")
        with open(cfg_path, "r", encoding="utf-8") as fh:
            override = yaml.safe_load(fh) or {}
        data = _deep_merge(data, override)

    return Config(data)
