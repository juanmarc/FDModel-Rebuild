"""Time-stepping methods for vorticity evolution."""

from __future__ import annotations

import numpy as np

from .fields import ModelState, state_from_vorticity
from .grid import PeriodicGrid
from .initial_conditions import zero_net_circulation
from .dynamics import vorticity_tendency


def rk4_step_vorticity(
    zeta: np.ndarray,
    grid: PeriodicGrid,
    dt: float,
    viscosity: float,
    *,
    jacobian_method: str = "arakawa",
    enforce_zero_circulation: bool = True,
) -> np.ndarray:
    """Advance vorticity by one fourth-order Runge-Kutta step."""

    if zeta.shape != grid.shape:
        raise ValueError(f"zeta shape {zeta.shape} does not match grid shape {grid.shape}")
    if dt <= 0.0:
        raise ValueError("dt must be positive")
    if viscosity < 0.0:
        raise ValueError("viscosity must be nonnegative")

    def rhs(field: np.ndarray) -> np.ndarray:
        return vorticity_tendency(field, grid, viscosity, jacobian_method=jacobian_method)

    k1 = rhs(zeta)
    k2 = rhs(zeta + 0.5 * dt * k1)
    k3 = rhs(zeta + 0.5 * dt * k2)
    k4 = rhs(zeta + dt * k3)

    next_zeta = zeta + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    if enforce_zero_circulation:
        next_zeta = zero_net_circulation(next_zeta)
    return next_zeta


def rk4_step_state(
    state: ModelState,
    dt: float,
    viscosity: float,
    *,
    jacobian_method: str = "arakawa",
    enforce_zero_circulation: bool = True,
) -> ModelState:
    """Advance a full model state by one RK4 vorticity step."""

    next_zeta = rk4_step_vorticity(
        state.zeta,
        state.grid,
        dt,
        viscosity,
        jacobian_method=jacobian_method,
        enforce_zero_circulation=enforce_zero_circulation,
    )
    return state_from_vorticity(next_zeta, state.grid)
