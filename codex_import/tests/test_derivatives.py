import numpy as np

from fdmodel.derivatives import d_dx, d_dy
from fdmodel.grid import make_periodic_grid


def test_centered_derivatives_are_periodic_across_boundary():
    grid = make_periodic_grid(64, 48)
    field = np.sin(grid.X) + 0.5 * np.cos(2.0 * grid.Y)

    expected_dx = np.cos(grid.X)
    expected_dy = -np.sin(2.0 * grid.Y)

    assert np.max(np.abs(d_dx(field, grid) - expected_dx)) < 3.0e-3
    assert np.max(np.abs(d_dy(field, grid) - expected_dy)) < 1.2e-2


def test_derivative_of_constant_is_zero_at_wrap():
    grid = make_periodic_grid(9, 7)
    field = np.full(grid.shape, 3.25)

    assert np.allclose(d_dx(field, grid), 0.0)
    assert np.allclose(d_dy(field, grid), 0.0)
