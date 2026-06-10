"""Plot plan-view and radial diagnostics from a saved state file."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from fdmodel.io import load_metadata_npz
from fdmodel.plotting import plot_state_file


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("state_file", help="NPZ state file created by runs/initialize_cases.py")
    parser.add_argument("--output-dir", default="runs/output/plots", help="directory for plot images")
    parser.add_argument("--bin-width", type=float, default=2.0e3, help="radial averaging bin width")
    parser.add_argument("--extent-radius", type=float, default=100.0e3, help="plan-view half-width")
    parser.add_argument("--max-radius", type=float, default=240.0e3, help="maximum radial-profile radius")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    os.environ.setdefault("MPLCONFIGDIR", str(PROJECT_ROOT / "runs" / "output" / "matplotlib-cache"))
    os.environ.setdefault("XDG_CACHE_HOME", str(PROJECT_ROOT / "runs" / "output" / "font-cache"))
    metadata = load_metadata_npz(args.state_file)
    title = f"{metadata['case']} {metadata['component']} initial conditions"
    plan_path, radial_path = plot_state_file(
        args.state_file,
        args.output_dir,
        title=title,
        bin_width=args.bin_width,
        extent_radius=args.extent_radius,
        max_radius=args.max_radius,
    )
    print(plan_path)
    print(radial_path)


if __name__ == "__main__":
    main()
