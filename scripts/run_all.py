from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def main(config_path: str) -> None:
    base = Path(__file__).resolve().parent
    scripts = [
        "exp1_overlap.py",
        "exp2_baseline_calibration.py",
        "exp3_prhs.py",
        "exp4_ge_vs_gecn.py",
    ]
    for script in scripts:
        cmd = [sys.executable, str(base / script), "--config", config_path]
        print("Running:", " ".join(cmd))
        subprocess.run(cmd, check=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    main(args.config)
