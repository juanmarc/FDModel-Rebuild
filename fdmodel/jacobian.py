"""Jacobian operators for the barotropic vorticity equation."""

from __future__ import annotations

import numpy as np

from .derivatives import d_dx, d_dy
from .grid import PeriodicGrid


def centered_jacobian(psi: np.ndarray, zeta: np.ndarray, grid: PeriodicGrid) -> np.ndarray:
    """Return ``J(psi, zeta) = psi_x * zeta_y - psi_y * zeta_x``.

    This is the straightforward second-order centered finite-difference form.
    """

    _check_shapes(psi, zeta, grid)
    return d_dx(psi, grid) * d_dy(zeta, grid) - d_dy(psi, grid) * d_dx(zeta, grid)


def arakawa_jacobian(psi: np.ndarray, zeta: np.ndarray, grid: PeriodicGrid) -> np.ndarray:
    """Return the second-order Arakawa Jacobian on a periodic grid.

    The Arakawa operator averages three centered stencil forms. It is preferred
    for long integrations because it has better discrete conservation
    properties than the plain centered Jacobian.
    """

    _check_shapes(psi, zeta, grid)
    dxdy = grid.dx * grid.dy

    psi_e = np.roll(psi, -1, axis=1)
    psi_w = np.roll(psi, 1, axis=1)
    psi_n = np.roll(psi, -1, axis=0)
    psi_s = np.roll(psi, 1, axis=0)
    psi_ne = np.roll(psi_e, -1, axis=0)
    psi_nw = np.roll(psi_w, -1, axis=0)
    psi_se = np.roll(psi_e, 1, axis=0)
    psi_sw = np.roll(psi_w, 1, axis=0)

    zeta_e = np.roll(zeta, -1, axis=1)
    zeta_w = np.roll(zeta, 1, axis=1)
    zeta_n = np.roll(zeta, -1, axis=0)
    zeta_s = np.roll(zeta, 1, axis=0)
    zeta_ne = np.roll(zeta_e, -1, axis=0)
    zeta_nw = np.roll(zeta_w, -1, axis=0)
    zeta_se = np.roll(zeta_e, 1, axis=0)
    zeta_sw = np.roll(zeta_w, 1, axis=0)

    j1 = ((psi_e - psi_w) * (zeta_n - zeta_s) - (psi_n - psi_s) * (zeta_e - zeta_w)) / (
        4.0 * dxdy
    )

    j2 = (
        psi_e * (zeta_ne - zeta_se)
        - psi_w * (zeta_nw - zeta_sw)
        - psi_n * (zeta_ne - zeta_nw)
        + psi_s * (zeta_se - zeta_sw)
    ) / (4.0 * dxdy)

    j3 = (
        zeta_n * (psi_ne - psi_nw)
        - zeta_s * (psi_se - psi_sw)
        - zeta_e * (psi_ne - psi_se)
        + zeta_w * (psi_nw - psi_sw)
    ) / (4.0 * dxdy)

    return (j1 + j2 + j3) / 3.0


def jacobian(
    psi: np.ndarray,
    zeta: np.ndarray,
    grid: PeriodicGrid,
    *,
    method: str = "arakawa",
) -> np.ndarray:
    """Return the requested Jacobian discretization."""

    if method == "arakawa":
        return arakawa_jacobian(psi, zeta, grid)
    if method == "centered":
        return centered_jacobian(psi, zeta, grid)
    raise ValueError(f"unknown Jacobian method {method!r}")


def _check_shapes(psi: np.ndarray, zeta: np.ndarray, grid: PeriodicGrid) -> None:
    if psi.shape != grid.shape:
        raise ValueError(f"psi shape {psi.shape} does not match grid shape {grid.shape}")
    if zeta.shape != grid.shape:
        raise ValueError(f"zeta shape {zeta.shape} does not match grid shape {grid.shape}")
