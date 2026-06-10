"""Generate initial-condition state files for monopole and ring cases."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from fdmodel.grid import make_centered_periodic_grid, make_periodic_grid
from fdmodel.initial_conditions import monopole_initial_states, schubert_ring_initial_states
from fdmodel.io import save_state_npz


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("monopole", "ring", "all"), help="case to initialize")
    parser.add_argument("--output-dir", default="runs/output", help="directory for NPZ state files")
    parser.add_argument("--nx", type=int, default=None, help="number of x grid points")
    parser.add_argument("--ny", type=int, default=None, help="number of y grid points")
    parser.add_argument("--dx", type=float, default=None, help="grid spacing in x")
    parser.add_argument("--dy", type=float, default=None, help="grid spacing in y")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    cases = ("monopole", "ring") if args.case == "all" else (args.case,)

    for case in cases:
        grid = grid_for_case(case, args.nx, args.ny, args.dx, args.dy)
        states = (
            monopole_initial_states(grid)
            if case == "monopole"
            else schubert_ring_initial_states(grid)
        )
        for component, state in states.items():
            path = output_dir / f"{case}_{component}.npz"
            save_state_npz(
                path,
                state,
                case=case,
                component=component,
                metadata={
                    "nx": grid.nx,
                    "ny": grid.ny,
                    "dx": grid.dx,
                    "dy": grid.dy,
                },
            )
            print(path)


def grid_for_case(
    case: str,
    nx: int | None,
    ny: int | None,
    dx: float | None,
    dy: float | None,
):
    if case == "monopole":
        default_nx = default_ny = 301
        default_dx = default_dy = 2.0e3
    else:
        default_nx = default_ny = 301
        default_dx = default_dy = 2.0e3

    nx = nx if nx is not None else default_nx
    ny = ny if ny is not None else default_ny
    dx = dx if dx is not None else default_dx
    dy = dy if dy is not None else default_dy
    if case == "monopole":
        return make_centered_periodic_grid(nx, ny, dx=dx, dy=dy)
    return make_periodic_grid(nx, ny, lx=nx * dx, ly=ny * dy)


if __name__ == "__main__":
    main()
