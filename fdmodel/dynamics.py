"""Right-hand-side tendencies for vorticity dynamics."""

from __future__ import annotations

import numpy as np

from .derivatives import laplacian
from .grid import PeriodicGrid
from .jacobian import jacobian
from .poisson import solve_poisson_fft


def vorticity_tendency(
    zeta: np.ndarray,
    grid: PeriodicGrid,
    viscosity: float,
    *,
    jacobian_method: str = "arakawa",
) -> np.ndarray:
    """Return ``d zeta / dt`` for the barotropic vorticity equation.

    The tendency is

    ``-J(psi, zeta) + viscosity * laplacian(zeta)``

    where ``psi`` is diagnosed from ``zeta`` with the periodic Poisson solver.
    """

    if zeta.shape != grid.shape:
        raise ValueError(f"zeta shape {zeta.shape} does not match grid shape {grid.shape}")
    if viscosity < 0.0:
        raise ValueError("viscosity must be nonnegative")

    psi = solve_poisson_fft(zeta, grid)
    advection = jacobian(psi, zeta, grid, method=jacobian_method)
    diffusion = laplacian(zeta, grid)
    return -advection + viscosity * diffusion
