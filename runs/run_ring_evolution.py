"""Run the hard-coded baseline ring-vortex evolution and save snapshots."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from fdmodel.evolution import EvolutionConfig, evolve_state
from fdmodel.grid import make_centered_periodic_grid
from fdmodel.initial_conditions import schubert_ring_initial_states
from fdmodel.io import save_state_npz


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="runs/output/ring_evolution")
    parser.add_argument("--max-time-count", type=int, default=86400)
    parser.add_argument("--output-interval", type=int, default=3600)
    parser.add_argument("--dt", type=float, default=1.0)
    parser.add_argument("--viscosity", type=float, default=100.0)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    grid = make_centered_periodic_grid(301, 301, dx=2.0e3, dy=2.0e3)
    initial_state = schubert_ring_initial_states(grid)["total"]
    config = EvolutionConfig(
        dt=args.dt,
        viscosity=args.viscosity,
        max_time_count=args.max_time_count,
        output_interval=args.output_interval,
    )
    max_speed = float((initial_state.u * initial_state.u + initial_state.v * initial_state.v).max() ** 0.5)
    advective_cfl = max_speed * config.dt / min(grid.dx, grid.dy)
    diffusion_number = config.viscosity * config.dt / (min(grid.dx, grid.dy) ** 2)
    if advective_cfl >= 1.0:
        raise ValueError(f"advective CFL must be less than 1; got {advective_cfl:.3f}")
    print(
        "run control: "
        f"dt={config.dt:g} s, dx={grid.dx:g} m, dy={grid.dy:g} m, "
        f"max_speed={max_speed:.6g} m/s, advective_cfl={advective_cfl:.6g}, "
        f"diffusion_number={diffusion_number:.6g}"
    )

    output_dir = Path(args.output_dir)
    for snapshot in evolve_state(initial_state, config):
        hour = snapshot.time_seconds / 3600.0
        path = output_dir / f"ring_step{snapshot.step:06d}_t{hour:06.2f}h.npz"
        save_state_npz(
            path,
            snapshot.state,
            case="ring",
            component="evolution",
            metadata={
                "step": snapshot.step,
                "time_seconds": snapshot.time_seconds,
                "time_hours": hour,
                "dt": config.dt,
                "viscosity": config.viscosity,
                "max_time_count": config.max_time_count,
                "output_interval": config.output_interval,
                "jacobian_method": config.jacobian_method,
                "initial_max_speed": max_speed,
                "initial_advective_cfl": advective_cfl,
                "diffusion_number": diffusion_number,
            },
        )
        print(path)


if __name__ == "__main__":
    main()
