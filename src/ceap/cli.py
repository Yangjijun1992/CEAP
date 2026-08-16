"""Command-line entry point for the CEAP afterpulse simulation."""
from __future__ import annotations

import argparse
import sys

from .config.loader import load_config
from .simulation.runner import run_simulation


def main(argv=None):
    parser = argparse.ArgumentParser(description="CEAP PMT afterpulse background simulation")
    parser.add_argument("--config", "-c", default="config/settings.yaml",
                       help="Simulation config YAML file (default: config/settings.yaml)")
    parser.add_argument("--run-id", help="Override run id")
    args = parser.parse_args(argv)

    cfg = load_config(args.config)
    if args.run_id:
        sim = cfg.get("simulation")
        sim["run_id"] = args.run_id

    result = run_simulation(cfg)
    print(f"run_id:            {result['run_id']}")
    print(f"n main pulses:     {result['n_main_pulses']} (muon: {result['n_muon_mains']})")
    print(f"windows scanned:   {result['n_windows_scanned']}")
    print(f"trigger events:    {result['n_trigger_events']}")
    print(f"background rate:   {result['background_rate_hz']:.4g} Hz")
    print(f"pdf:               {result['pdf']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
