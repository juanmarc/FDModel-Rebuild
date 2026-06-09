import numpy as np

from fdmodel.grid import make_periodic_grid


def test_periodic_grid_shape_and_spacing():
    grid = make_periodic_grid(8, 6, lx=4.0, ly=3.0)

    assert grid.shape == (6, 8)
    assert grid.X.shape == (6, 8)
    assert grid.Y.shape == (6, 8)
    assert np.isclose(grid.dx, 0.5)
    assert np.isclose(grid.dy, 0.5)
    assert np.isclose(grid.x[-1], 3.5)
    assert np.isclose(grid.y[-1], 2.5)
