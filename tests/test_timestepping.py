import numpy as np
import pytest

from fdmodel.fields import state_from_vorticity
from fdmodel.grid import make_periodic_grid
from fdmodel.timestepping import rk4_step_state, rk4_step_vorticity


def test_rk4_zero_vorticity_remains_zero():
    grid = make_periodic_grid(32, 32)
    zeta = np.zeros(grid.shape)

    next_zeta = rk4_step_vorticity(zeta, grid, dt=0.1, viscosity=0.01)

    assert np.allclose(next_zeta, 0.0)


def test_rk4_step_preserves_shape_and_zero_mean():
    grid = make_periodic_grid(48, 44)
    zeta = np.sin(grid.X) + 0.2 * np.cos(2.0 * grid.Y)
    zeta = zeta - np.mean(zeta)

    next_zeta = rk4_step_vorticity(zeta, grid, dt=0.01, viscosity=0.001)

    assert next_zeta.shape == grid.shape
    assert abs(np.mean(next_zeta)) < 1.0e-15
    assert np.all(np.isfinite(next_zeta))


def test_rk4_diffusion_of_single_fourier_mode_matches_short_exact_solution():
    grid = make_periodic_grid(96, 96)
    zeta = np.sin(2.0 * grid.X) * np.cos(3.0 * grid.Y)
    viscosity = 0.01
    dt = 0.02
    wavenumber_squared = 2.0**2 + 3.0**2

    next_zeta = rk4_step_vorticity(zeta, grid, dt=dt, viscosity=viscosity)
    expected = np.exp(-viscosity * wavenumber_squared * dt) * zeta

    assert np.max(np.abs(next_zeta - expected)) < 2.0e-5


def test_rk4_step_state_recomputes_diagnostic_fields():
    grid = make_periodic_grid(40, 36)
    zeta = np.sin(grid.X) * np.cos(grid.Y)
    state = state_from_vorticity(zeta, grid)

    next_state = rk4_step_state(state, dt=0.01, viscosity=0.001)

    assert next_state.grid is grid
    assert next_state.zeta.shape == grid.shape
    assert next_state.psi.shape == grid.shape
    assert next_state.u.shape == grid.shape
    assert next_state.v.shape == grid.shape


def test_rk4_rejects_invalid_inputs():
    grid = make_periodic_grid(20, 18)
    zeta = np.zeros(grid.shape)

    with pytest.raises(ValueError, match="dt must be positive"):
        rk4_step_vorticity(zeta, grid, dt=0.0, viscosity=0.0)
    with pytest.raises(ValueError, match="viscosity"):
        rk4_step_vorticity(zeta, grid, dt=1.0, viscosity=-1.0)
    with pytest.raises(ValueError, match="zeta shape"):
        rk4_step_vorticity(np.zeros((grid.ny + 1, grid.nx)), grid, dt=1.0, viscosity=0.0)
