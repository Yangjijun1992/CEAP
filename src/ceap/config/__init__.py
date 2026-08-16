"""Configuration loading and management.

Detector parameters are reserved for later provisioning via YAML config files.
"""
from .loader import load_config, DEFAULTS_DIR

__all__ = ["load_config", "DEFAULTS_DIR"]
