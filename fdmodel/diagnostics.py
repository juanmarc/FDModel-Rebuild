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


def azimuthal_mean_profile(
    field: np.ndarray,
    grid: PeriodicGrid,
    bin_width: float,
    center: tuple[float, float] | None = None,
    max_radius: float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return radial bin centers, azimuthal means, and sample counts."""

    if field.shape != grid.shape:
        raise ValueError(f"field shape {field.shape} does not match grid shape {grid.shape}")
    if bin_width <= 0.0:
        raise ValueError("bin_width must be positive")

    _, _, radius = _radial_distance(grid, center)
    if max_radius is None:
        max_radius = float(np.max(radius))
    nbins = int(np.floor(max_radius / bin_width)) + 1
    bin_index = np.floor(radius / bin_width).astype(int)
    mask = bin_index < nbins

    sums = np.bincount(bin_index[mask].ravel(), weights=field[mask].ravel(), minlength=nbins)
    counts = np.bincount(bin_index[mask].ravel(), minlength=nbins)
    means = np.full(nbins, np.nan, dtype=float)
    valid = counts > 0
    means[valid] = sums[valid] / counts[valid]
    centers = (np.arange(nbins) + 0.5) * bin_width
    return centers, means, counts


def radial_wind_profile(
    u: np.ndarray,
    v: np.ndarray,
    grid: PeriodicGrid,
    bin_width: float,
    center: tuple[float, float] | None = None,
    max_radius: float | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return azimuthal-mean tangential wind as a radial profile."""

    if u.shape != grid.shape or v.shape != grid.shape:
        raise ValueError("u and v must both match grid shape")

    dx, dy, radius = _radial_distance(grid, center)
    tangential = np.zeros(grid.shape, dtype=float)
    nonzero = radius > 0.0
    tangential[nonzero] = (-u[nonzero] * dy[nonzero] + v[nonzero] * dx[nonzero]) / radius[nonzero]
    return azimuthal_mean_profile(tangential, grid, bin_width, center=center, max_radius=max_radius)


def _radial_distance(
    grid: PeriodicGrid,
    center: tuple[float, float] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    cx, cy = center if center is not None else (0.5 * grid.lx, 0.5 * grid.ly)
    dx = (grid.X - cx + 0.5 * grid.lx) % grid.lx - 0.5 * grid.lx
    dy = (grid.Y - cy + 0.5 * grid.ly) % grid.ly - 0.5 * grid.ly
    return dx, dy, np.hypot(dx, dy)
