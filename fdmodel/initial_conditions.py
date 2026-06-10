"""Initial vorticity fields for periodic domains."""

from __future__ import annotations

import numpy as np

from .fields import ModelState, state_from_vorticity
from .grid import PeriodicGrid


def periodic_displacement(coord: np.ndarray, center: float, length: float) -> np.ndarray:
    """Return minimum-image displacement from ``center`` on a periodic interval."""

    return (coord - center + 0.5 * length) % length - 0.5 * length


def zero_net_circulation(zeta: np.ndarray, grid: PeriodicGrid | None = None) -> np.ndarray:
    """Subtract a constant so area-integrated circulation is zero."""

    if grid is None or not grid.endpoint_inclusive:
        correction = float(np.mean(zeta))
    else:
        if zeta.shape != grid.shape:
            raise ValueError(f"zeta shape {zeta.shape} does not match grid shape {grid.shape}")
        area = grid.lx * grid.ly
        integral = float(np.trapezoid(np.trapezoid(zeta, x=grid.x, axis=1), x=grid.y, axis=0))
        correction = integral / area
    return zeta - correction


def radial_distance(
    grid: PeriodicGrid,
    center: tuple[float, float] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return minimum-image x/y displacement and radius from ``center``."""

    cx, cy = center if center is not None else grid.center
    dx = periodic_displacement(grid.X, cx, grid.lx)
    dy = periodic_displacement(grid.Y, cy, grid.ly)
    return dx, dy, np.hypot(dx, dy)


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

    cx, cy = center if center is not None else grid.center
    dx = periodic_displacement(grid.X, cx, grid.lx)
    dy = periodic_displacement(grid.Y, cy, grid.ly)
    zeta = amplitude * np.exp(-(dx * dx + dy * dy) / (2.0 * sigma * sigma))

    if enforce_zero_circulation:
        zeta = zero_net_circulation(zeta, grid)
    return zeta


def monopole_base_vorticity(
    grid: PeriodicGrid,
    zeta_max: float = 1.0e-3,
    decay_radius: float = 47.0e3,
    center: tuple[float, float] | None = None,
    enforce_zero_circulation: bool = True,
) -> np.ndarray:
    """Return the Gaussian monopole basic-state vorticity.

    The profile is ``zeta_max * exp(-(r / decay_radius)**2)``. Distances are in
    the same units as the grid, normally meters for physical runs.
    """

    if zeta_max <= 0.0:
        raise ValueError("zeta_max must be positive")
    if decay_radius <= 0.0:
        raise ValueError("decay_radius must be positive")

    _, _, radius = radial_distance(grid, center)
    zeta = zeta_max * np.exp(-((radius / decay_radius) ** 2))

    if enforce_zero_circulation:
        zeta = zero_net_circulation(zeta, grid)
    return zeta


def monopole_perturbation_vorticity(
    grid: PeriodicGrid,
    zeta_max: float = 0.5e-3,
    decay_radius: float = (2.0 / 3.0) * 47.0e3,
    center: tuple[float, float] | None = None,
    enforce_zero_circulation: bool = True,
) -> np.ndarray:
    """Return a localized Gaussian monopole perturbation."""

    if center is None:
        cx, cy = grid.center
        center = (cx + 54.0e3, cy)
    perturbation = monopole_base_vorticity(
        grid,
        zeta_max=zeta_max,
        decay_radius=decay_radius,
        center=center,
        enforce_zero_circulation=False,
    )

    if enforce_zero_circulation:
        perturbation = zero_net_circulation(perturbation, grid)
    return perturbation


def monopole_initial_states(
    grid: PeriodicGrid,
    zeta_max: float = 1.0e-3,
    decay_radius: float = 47.0e3,
    perturbation_fraction: float = 0.5,
    perturbation_decay_fraction: float = 2.0 / 3.0,
    perturbation_center: tuple[float, float] | None = None,
) -> dict[str, ModelState]:
    """Initialize base, perturbation, and total monopole states."""

    base = monopole_base_vorticity(grid, zeta_max=zeta_max, decay_radius=decay_radius)
    perturbation = monopole_perturbation_vorticity(
        grid,
        zeta_max=perturbation_fraction * zeta_max,
        decay_radius=perturbation_decay_fraction * decay_radius,
        center=perturbation_center,
    )
    total = zero_net_circulation(base + perturbation, grid)
    return {
        "base": state_from_vorticity(base, grid),
        "perturbation": state_from_vorticity(perturbation, grid),
        "total": state_from_vorticity(total, grid),
    }


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

    if center is None:
        grid_cx, grid_cy = grid.center
        center = (grid_cx + 0.1 * grid.lx, grid_cy)
    cx, cy = center
    dx = periodic_displacement(grid.X, cx, grid.lx)
    dy = periodic_displacement(grid.Y, cy, grid.ly)
    perturbation = amplitude * np.exp(-(dx * dx + dy * dy) / (2.0 * sigma * sigma))
    result = zeta + perturbation

    if enforce_zero_circulation:
        result = zero_net_circulation(result, grid)
    return result


def schubert_ring_base_vorticity(
    grid: PeriodicGrid,
    inner_core_radius: float = 15.0e3,
    inner_ramp_radius: float = 22.5e3,
    outer_ramp_start_radius: float = 25.0e3,
    outer_radius: float = 32.5e3,
    core_vorticity: float = 4.1825e-4,
    ring_vorticity: float = 7.000808e-3,
    center: tuple[float, float] | None = None,
    enforce_zero_circulation: bool = True,
) -> np.ndarray:
    """Return the axisymmetric Schubert-style ring-vortex profile.

    This follows the piecewise cubic profile in ``FDModel-main/basic.state.f``
    with the azimuthal perturbation omitted.
    """

    _, _, radius = radial_distance(grid, center)
    zeta = np.zeros(grid.shape, dtype=float)

    zeta[radius <= inner_core_radius] = core_vorticity

    mask = (radius > inner_core_radius) & (radius <= inner_ramp_radius)
    fn1 = (radius[mask] - inner_core_radius) / (inner_ramp_radius - inner_core_radius)
    fn2 = (inner_ramp_radius - radius[mask]) / (inner_ramp_radius - inner_core_radius)
    s1 = smooth_step_down(fn1)
    s2 = smooth_step_down(fn2)
    zeta[mask] = core_vorticity * s1 + ring_vorticity * s2

    mask = (radius > inner_ramp_radius) & (radius <= outer_ramp_start_radius)
    zeta[mask] = ring_vorticity

    mask = (radius > outer_ramp_start_radius) & (radius <= outer_radius)
    fn = (radius[mask] - outer_ramp_start_radius) / (outer_radius - outer_ramp_start_radius)
    zeta[mask] = ring_vorticity * smooth_step_down(fn)

    if enforce_zero_circulation:
        zeta = zero_net_circulation(zeta, grid)
    return zeta


def schubert_ring_perturbation_vorticity(
    grid: PeriodicGrid,
    perturbation_amplitude: float = 1.0e-5,
    inner_ramp_radius: float = 22.5e3,
    outer_ramp_start_radius: float = 25.0e3,
    outer_radius: float = 32.5e3,
    max_wavenumber: int = 8,
    center: tuple[float, float] | None = None,
    enforce_zero_circulation: bool = True,
) -> np.ndarray:
    """Return the azimuthal perturbation used by the legacy ring initializer."""

    if max_wavenumber < 1:
        raise ValueError("max_wavenumber must be positive")

    dx, dy, radius = radial_distance(grid, center)
    theta = np.mod(np.arctan2(dy, dx), 2.0 * np.pi)
    sum_cos = np.zeros(grid.shape, dtype=float)
    for wavenumber in range(1, max_wavenumber + 1):
        sum_cos += np.cos(wavenumber * theta)

    envelope = np.zeros(grid.shape, dtype=float)

    mask = (radius > 15.0e3) & (radius <= inner_ramp_radius)
    fn = (inner_ramp_radius - radius[mask] + 1.0) / (inner_ramp_radius - 15.0e3)
    envelope[mask] = smooth_step_down(fn)

    mask = (radius > inner_ramp_radius) & (radius <= outer_ramp_start_radius)
    envelope[mask] = 1.0

    mask = (radius > outer_ramp_start_radius) & (radius <= outer_radius)
    fn = (radius[mask] - outer_ramp_start_radius) / (outer_radius - outer_ramp_start_radius)
    envelope[mask] = smooth_step_down(fn)

    perturbation = perturbation_amplitude * sum_cos * envelope
    if enforce_zero_circulation:
        perturbation = zero_net_circulation(perturbation, grid)
    return perturbation


def schubert_ring_initial_states(grid: PeriodicGrid) -> dict[str, ModelState]:
    """Initialize base, perturbation, and total ring-vortex states."""

    base = schubert_ring_base_vorticity(grid)
    perturbation = schubert_ring_perturbation_vorticity(grid)
    total = zero_net_circulation(base + perturbation, grid)
    return {
        "base": state_from_vorticity(base, grid),
        "perturbation": state_from_vorticity(perturbation, grid),
        "total": state_from_vorticity(total, grid),
    }


def smooth_step_down(fraction: np.ndarray) -> np.ndarray:
    """Return ``1 - 3 f**2 + 2 f**3`` for cubic radial transitions."""

    return 1.0 - 3.0 * fraction * fraction + 2.0 * fraction * fraction * fraction
