"""FFT-based Poisson inversion on periodic grids."""

from __future__ import annotations

import numpy as np

from .grid import PeriodicGrid


def solve_poisson_fft(zeta: np.ndarray, grid: PeriodicGrid) -> np.ndarray:
    """Solve ``laplacian(psi) = zeta`` with periodic boundaries.

    The zero Fourier mode is set to zero, fixing the arbitrary additive
    constant in ``psi``. The input vorticity is mean-corrected before inversion
    because periodic Poisson problems require a zero-mean right-hand side.
    """

    if zeta.shape != grid.shape:
        raise ValueError(f"zeta shape {zeta.shape} does not match grid shape {grid.shape}")

    rhs = zeta - np.mean(zeta)
    rhs_hat = np.fft.fft2(rhs)

    kx = 2.0 * np.pi * np.fft.fftfreq(grid.nx, d=grid.dx)
    ky = 2.0 * np.pi * np.fft.fftfreq(grid.ny, d=grid.dy)
    kx2, ky2 = np.meshgrid(kx * kx, ky * ky, indexing="xy")
    k2 = kx2 + ky2

    psi_hat = np.zeros_like(rhs_hat, dtype=np.complex128)
    nonzero = k2 > 0.0
    psi_hat[nonzero] = -rhs_hat[nonzero] / k2[nonzero]
    return np.fft.ifft2(psi_hat).real
