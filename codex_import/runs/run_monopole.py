"""Build and diagnose a Gaussian monopole initial condition."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from fdmodel.diagnostics import summarize_vorticity, velocity_from_streamfunction
from fdmodel.grid import make_periodic_grid
from fdmodel.initial_conditions import add_gaussian_perturbation, gaussian_monopole_vorticity
from fdmodel.poisson import solve_poisson_fft


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nx", type=int, default=128, help="number of grid points in x")
    parser.add_argument("--ny", type=int, default=128, help="number of grid points in y")
    parser.add_argument("--lx", type=float, default=6.283185307179586, help="domain length in x")
    parser.add_argument("--ly", type=float, default=6.283185307179586, help="domain length in y")
    parser.add_argument("--amplitude", type=float, default=1.0, help="monopole vorticity amplitude")
    parser.add_argument("--sigma", type=float, default=0.4, help="monopole Gaussian width")
    parser.add_argument("--perturbation-amplitude", type=float, default=0.05, help="perturbation amplitude")
    parser.add_argument("--perturbation-sigma", type=float, default=0.2, help="perturbation Gaussian width")
    parser.add_argument("--no-perturbation", action="store_true", help="skip localized perturbation")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    grid = make_periodic_grid(args.nx, args.ny, args.lx, args.ly)

    zeta = gaussian_monopole_vorticity(grid, amplitude=args.amplitude, sigma=args.sigma)
    if not args.no_perturbation:
        zeta = add_gaussian_perturbation(
            zeta,
            grid,
            amplitude=args.perturbation_amplitude,
            sigma=args.perturbation_sigma,
        )

    psi = solve_poisson_fft(zeta, grid)
    u, v = velocity_from_streamfunction(psi, grid)

    diagnostics = summarize_vorticity(zeta, grid)
    diagnostics["max_speed"] = float((u * u + v * v).max() ** 0.5)

    print(f"grid: nx={grid.nx}, ny={grid.ny}, dx={grid.dx:.6g}, dy={grid.dy:.6g}")
    for name, value in diagnostics.items():
        print(f"{name}: {value:.12g}")


if __name__ == "__main__":
    main()
