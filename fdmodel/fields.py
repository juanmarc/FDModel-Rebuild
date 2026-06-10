"""Containers and constructors for model state fields."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .diagnostics import velocity_from_streamfunction
from .grid import PeriodicGrid
from .poisson import solve_poisson_fft


@dataclass(frozen=True)
class ModelState:
    """Vorticity, streamfunction, and velocity on a common grid."""

    grid: PeriodicGrid
    zeta: np.ndarray
    psi: np.ndarray
    u: np.ndarray
    v: np.ndarray


def state_from_vorticity(zeta: np.ndarray, grid: PeriodicGrid) -> ModelState:
    """Invert vorticity and compute velocity diagnostics."""

    if zeta.shape != grid.shape:
        raise ValueError(f"zeta shape {zeta.shape} does not match grid shape {grid.shape}")

    psi = solve_poisson_fft(zeta, grid)
    u, v = velocity_from_streamfunction(psi, grid)
    return ModelState(grid=grid, zeta=zeta, psi=psi, u=u, v=v)
