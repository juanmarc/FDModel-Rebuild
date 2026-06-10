import numpy as np
import pytest

from fdmodel.grid import make_periodic_grid
from fdmodel.jacobian import arakawa_jacobian, centered_jacobian, jacobian


def test_centered_jacobian_matches_analytic_periodic_fields():
    grid = make_periodic_grid(160, 144)
    psi = np.sin(2.0 * grid.X) * np.cos(grid.Y)
    zeta = np.cos(grid.X) * np.sin(3.0 * grid.Y)

    expected = (
        6.0 * np.cos(2.0 * grid.X) * np.cos(grid.Y) * np.cos(grid.X) * np.cos(3.0 * grid.Y)
        - np.sin(2.0 * grid.X) * np.sin(grid.Y) * np.sin(grid.X) * np.sin(3.0 * grid.Y)
    )

    actual = centered_jacobian(psi, zeta, grid)

    assert np.max(np.abs(actual - expected)) < 2.5e-2


def test_arakawa_jacobian_matches_analytic_periodic_fields():
    grid = make_periodic_grid(160, 144)
    psi = np.sin(2.0 * grid.X) * np.cos(grid.Y)
    zeta = np.cos(grid.X) * np.sin(3.0 * grid.Y)

    expected = (
        6.0 * np.cos(2.0 * grid.X) * np.cos(grid.Y) * np.cos(grid.X) * np.cos(3.0 * grid.Y)
        - np.sin(2.0 * grid.X) * np.sin(grid.Y) * np.sin(grid.X) * np.sin(3.0 * grid.Y)
    )

    actual = arakawa_jacobian(psi, zeta, grid)

    assert np.max(np.abs(actual - expected)) < 3.5e-2


def test_jacobian_of_field_with_itself_is_zero():
    grid = make_periodic_grid(49, 45)
    field = np.sin(grid.X) + 0.3 * np.cos(2.0 * grid.Y)

    assert np.allclose(centered_jacobian(field, field, grid), 0.0)
    assert np.allclose(arakawa_jacobian(field, field, grid), 0.0)


def test_arakawa_jacobian_has_near_zero_domain_integral():
    grid = make_periodic_grid(63, 57)
    psi = np.sin(grid.X) * np.cos(2.0 * grid.Y)
    zeta = np.cos(3.0 * grid.X) + 0.5 * np.sin(grid.Y)

    domain_integral = np.sum(arakawa_jacobian(psi, zeta, grid)) * grid.cell_area

    assert abs(domain_integral) < 1.0e-13


def test_jacobian_dispatch_and_invalid_method():
    grid = make_periodic_grid(16, 12)
    psi = np.sin(grid.X)
    zeta = np.cos(grid.Y)

    assert np.allclose(jacobian(psi, zeta, grid), arakawa_jacobian(psi, zeta, grid))
    assert np.allclose(jacobian(psi, zeta, grid, method="centered"), centered_jacobian(psi, zeta, grid))
    with pytest.raises(ValueError, match="unknown Jacobian method"):
        jacobian(psi, zeta, grid, method="bad")


def test_jacobian_rejects_shape_mismatch():
    grid = make_periodic_grid(16, 12)
    psi = np.zeros(grid.shape)
    zeta = np.zeros((grid.ny + 1, grid.nx))

    with pytest.raises(ValueError, match="zeta shape"):
        arakawa_jacobian(psi, zeta, grid)
