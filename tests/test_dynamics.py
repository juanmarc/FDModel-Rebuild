import numpy as np
import pytest

from fdmodel.derivatives import laplacian
from fdmodel.dynamics import vorticity_tendency
from fdmodel.grid import make_periodic_grid
from fdmodel.jacobian import jacobian
from fdmodel.poisson import solve_poisson_fft


def test_zero_vorticity_has_zero_tendency():
    grid = make_periodic_grid(32, 24)
    zeta = np.zeros(grid.shape)

    tendency = vorticity_tendency(zeta, grid, viscosity=100.0)

    assert np.allclose(tendency, 0.0)


def test_single_fourier_mode_reduces_to_diffusion():
    grid = make_periodic_grid(96, 96)
    zeta = np.sin(2.0 * grid.X) * np.cos(3.0 * grid.Y)
    viscosity = 0.01

    tendency = vorticity_tendency(zeta, grid, viscosity=viscosity)
    expected = viscosity * laplacian(zeta, grid)

    assert np.max(np.abs(tendency - expected)) < 1.0e-13


def test_inviscid_tendency_has_near_zero_domain_mean():
    grid = make_periodic_grid(64, 60)
    zeta = np.sin(grid.X) + 0.25 * np.cos(2.0 * grid.X + grid.Y)

    tendency = vorticity_tendency(zeta, grid, viscosity=0.0)

    assert abs(np.mean(tendency)) < 1.0e-15


def test_tendency_accepts_centered_jacobian_method():
    grid = make_periodic_grid(40, 36)
    zeta = np.sin(grid.X) + 0.1 * np.cos(grid.Y)

    tendency = vorticity_tendency(zeta, grid, viscosity=0.0, jacobian_method="centered")

    assert tendency.shape == grid.shape
    assert np.all(np.isfinite(tendency))


def test_tendency_rejects_invalid_inputs():
    grid = make_periodic_grid(20, 18)
    zeta = np.zeros((grid.ny + 1, grid.nx))

    with pytest.raises(ValueError, match="zeta shape"):
        vorticity_tendency(zeta, grid, viscosity=0.0)
    with pytest.raises(ValueError, match="viscosity"):
        vorticity_tendency(np.zeros(grid.shape), grid, viscosity=-1.0)


def test_tendency_matches_manual_rhs_composition():
    grid = make_periodic_grid(48, 44)
    zeta = np.sin(grid.X) + 0.2 * np.cos(2.0 * grid.Y)
    viscosity = 0.05

    psi = solve_poisson_fft(zeta, grid)
    tendency = vorticity_tendency(zeta, grid, viscosity=viscosity)
    expected = -jacobian(psi, zeta, grid) + viscosity * laplacian(zeta, grid)

    assert psi.shape == grid.shape
    assert np.max(np.abs(tendency - expected)) < 1.0e-13
