"""Initial vorticity fields for periodic domains."""

from __future__ import annotations

import numpy as np

from .grid import PeriodicGrid


def periodic_displacement(coord: np.ndarray, center: float, length: float) -> np.ndarray:
    """Return minimum-image displacement from ``center`` on a periodic interval."""

    return (coord - center + 0.5 * length) % length - 0.5 * length


def zero_net_circulation(zeta: np.ndarray) -> np.ndarray:
    """Subtract the spatial mean so total circulation is zero on a uniform grid."""

    return zeta - np.mean(zeta)


def gaussian_monopole_vorticity(
    grid: PeriodicGrid,
    amplitude: float = 1.0,
    sigma: float = 0.4,
    center: tuple[float, float] | None = None,
    enforce_zero_circulation: bool = True,
) -> np.ndarray:
    """Create a Gaussian monopole vorticity field on a periodic grid."""

    if sigma <= 0.0:
        raise ValueError("sigma must be positive")

    cx, cy = center if center is not None else (0.5 * grid.lx, 0.5 * grid.ly)
    dx = periodic_displacement(grid.X, cx, grid.lx)
    dy = periodic_displacement(grid.Y, cy, grid.ly)
    zeta = amplitude * np.exp(-(dx * dx + dy * dy) / (2.0 * sigma * sigma))

    if enforce_zero_circulation:
        zeta = zero_net_circulation(zeta)
    return zeta


def add_gaussian_perturbation(
    zeta: np.ndarray,
    grid: PeriodicGrid,
    amplitude: float = 0.1,
    sigma: float = 0.2,
    center: tuple[float, float] | None = None,
    enforce_zero_circulation: bool = True,
) -> np.ndarray:
    """Add a localized Gaussian perturbation to an existing vorticity field."""

    if zeta.shape != grid.shape:
        raise ValueError(f"zeta shape {zeta.shape} does not match grid shape {grid.shape}")
    if sigma <= 0.0:
        raise ValueError("sigma must be positive")

    cx, cy = center if center is not None else (0.6 * grid.lx, 0.5 * grid.ly)
    dx = periodic_displacement(grid.X, cx, grid.lx)
    dy = periodic_displacement(grid.Y, cy, grid.ly)
    perturbation = amplitude * np.exp(-(dx * dx + dy * dy) / (2.0 * sigma * sigma))
    result = zeta + perturbation

    if enforce_zero_circulation:
        result = zero_net_circulation(result)
    return result
