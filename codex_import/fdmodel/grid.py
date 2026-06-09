"""Periodic Cartesian grids."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class PeriodicGrid:
    """Uniform two-dimensional periodic Cartesian grid.

    The coordinate arrays use ``endpoint=False`` so the last grid point is one
    spacing before the periodic wrap back to zero.
    """

    nx: int
    ny: int
    lx: float
    ly: float
    x: np.ndarray
    y: np.ndarray
    X: np.ndarray
    Y: np.ndarray

    @property
    def dx(self) -> float:
        return self.lx / self.nx

    @property
    def dy(self) -> float:
        return self.ly / self.ny

    @property
    def shape(self) -> tuple[int, int]:
        return (self.ny, self.nx)

    @property
    def cell_area(self) -> float:
        return self.dx * self.dy


def make_periodic_grid(
    nx: int,
    ny: int,
    lx: float = 2.0 * np.pi,
    ly: float = 2.0 * np.pi,
) -> PeriodicGrid:
    """Create a uniform periodic grid with arrays shaped ``(ny, nx)``."""

    if nx <= 0 or ny <= 0:
        raise ValueError("nx and ny must be positive")
    if lx <= 0.0 or ly <= 0.0:
        raise ValueError("domain lengths must be positive")

    x = np.linspace(0.0, lx, nx, endpoint=False)
    y = np.linspace(0.0, ly, ny, endpoint=False)
    X, Y = np.meshgrid(x, y, indexing="xy")
    return PeriodicGrid(nx=nx, ny=ny, lx=float(lx), ly=float(ly), x=x, y=y, X=X, Y=Y)
