"""Basic scalar diagnostics for vorticity and streamfunction fields."""

from __future__ import annotations

import numpy as np

from .derivatives import d_dx, d_dy
from .grid import PeriodicGrid


def velocity_from_streamfunction(psi: np.ndarray, grid: PeriodicGrid) -> tuple[np.ndarray, np.ndarray]:
    """Return velocity components ``u = -dpsi/dy`` and ``v = dpsi/dx``."""

    return -d_dy(psi, grid), d_dx(psi, grid)


def circulation(zeta: np.ndarray, grid: PeriodicGrid) -> float:
    """Return total circulation, the area integral of vorticity."""

    if zeta.shape != grid.shape:
        raise ValueError(f"zeta shape {zeta.shape} does not match grid shape {grid.shape}")
    return float(np.sum(zeta) * grid.cell_area)


def enstrophy(zeta: np.ndarray, grid: PeriodicGrid) -> float:
    """Return enstrophy, ``0.5 * integral(zeta**2 dA)``."""

    if zeta.shape != grid.shape:
        raise ValueError(f"zeta shape {zeta.shape} does not match grid shape {grid.shape}")
    return float(0.5 * np.sum(zeta * zeta) * grid.cell_area)


def vorticity_extrema(zeta: np.ndarray) -> tuple[float, float]:
    """Return minimum and maximum vorticity."""

    return float(np.min(zeta)), float(np.max(zeta))


def summarize_vorticity(zeta: np.ndarray, grid: PeriodicGrid) -> dict[str, float]:
    """Return a compact diagnostics dictionary for a vorticity field."""

    zeta_min, zeta_max = vorticity_extrema(zeta)
    return {
        "circulation": circulation(zeta, grid),
        "zeta_min": zeta_min,
        "zeta_max": zeta_max,
        "enstrophy": enstrophy(zeta, grid),
    }
