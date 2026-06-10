import numpy as np

import pytest

from fdmodel.grid import make_centered_periodic_grid, make_periodic_grid


def test_periodic_grid_shape_and_spacing():
    grid = make_periodic_grid(8, 6, lx=4.0, ly=3.0)

    assert grid.shape == (6, 8)
    assert grid.X.shape == (6, 8)
    assert grid.Y.shape == (6, 8)
    assert np.isclose(grid.dx, 0.5)
    assert np.isclose(grid.dy, 0.5)
    assert np.isclose(grid.x[-1], 3.5)
    assert np.isclose(grid.y[-1], 2.5)


def test_centered_periodic_grid_has_true_center_point():
    grid = make_centered_periodic_grid(301, 301, dx=2.0e3, dy=2.0e3)

    assert grid.shape == (301, 301)
    assert grid.center == (0.0, 0.0)
    assert np.isclose(grid.lx, 600.0e3)
    assert np.isclose(grid.ly, 600.0e3)
    assert np.isclose(grid.dx, 2.0e3)
    assert np.isclose(grid.dy, 2.0e3)
    assert np.isclose(grid.x[150], 0.0)
    assert np.isclose(grid.y[150], 0.0)
    assert np.isclose(grid.x[0], -300.0e3)
    assert np.isclose(grid.x[-1], 300.0e3)


def test_centered_periodic_grid_requires_odd_dimensions():
    with pytest.raises(ValueError, match="odd"):
        make_centered_periodic_grid(300, 301, dx=2.0e3, dy=2.0e3)
