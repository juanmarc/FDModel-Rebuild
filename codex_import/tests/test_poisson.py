import numpy as np

from fdmodel.grid import make_periodic_grid
from fdmodel.poisson import solve_poisson_fft


def test_poisson_fft_recovers_known_streamfunction_up_to_constant():
    grid = make_periodic_grid(64, 64)
    psi_exact = np.sin(2.0 * grid.X) * np.cos(3.0 * grid.Y)
    zeta = -(2.0**2 + 3.0**2) * psi_exact

    psi = solve_poisson_fft(zeta, grid)

    assert abs(np.mean(psi)) < 1.0e-14
    assert np.max(np.abs(psi - psi_exact)) < 1.0e-12


def test_poisson_fft_handles_nonzero_mean_rhs_by_removing_zero_mode():
    grid = make_periodic_grid(32, 32)
    zeta = -np.sin(grid.X) + 7.0

    psi = solve_poisson_fft(zeta, grid)

    assert abs(np.mean(psi)) < 1.0e-14
    assert np.max(np.abs(psi - np.sin(grid.X))) < 1.0e-12
