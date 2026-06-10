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
    center_x: float | None = None
    center_y: float | None = None

    @property
    def dx(self) -> float:
        return float(abs(self.x[1] - self.x[0])) if self.nx > 1 else self.lx

    @property
    def dy(self) -> float:
        return float(abs(self.y[1] - self.y[0])) if self.ny > 1 else self.ly

    @property
    def shape(self) -> tuple[int, int]:
        return (self.ny, self.nx)

    @property
    def cell_area(self) -> float:
        return self.dx * self.dy

    @property
    def endpoint_inclusive(self) -> bool:
        """Return True when coordinates include both physical domain endpoints."""

        x_span = abs(float(self.x[-1] - self.x[0])) if self.nx > 1 else 0.0
        y_span = abs(float(self.y[-1] - self.y[0])) if self.ny > 1 else 0.0
        return np.isclose(x_span, self.lx) and np.isclose(y_span, self.ly)

    @property
    def center(self) -> tuple[float, float]:
        """Return the physical center used for radial diagnostics."""

        cx = self.center_x if self.center_x is not None else 0.5 * self.lx
        cy = self.center_y if self.center_y is not None else 0.5 * self.ly
        return (float(cx), float(cy))


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


def make_centered_periodic_grid(
    nx: int,
    ny: int,
    dx: float,
    dy: float,
) -> PeriodicGrid:
    """Create an odd endpoint-inclusive grid with a true center point.

    This follows the legacy Fortran convention where ``nx`` grid points span
    ``nx - 1`` intervals. For example, 301 points at 2 km spacing represent a
    600 km domain with coordinates from -300 km through +300 km.
    """

    if nx <= 0 or ny <= 0:
        raise ValueError("nx and ny must be positive")
    if nx % 2 == 0 or ny % 2 == 0:
        raise ValueError("centered grids require odd nx and ny")
    if dx <= 0.0 or dy <= 0.0:
        raise ValueError("grid spacing must be positive")

    x = (np.arange(nx, dtype=float) - nx // 2) * dx
    y = (np.arange(ny, dtype=float) - ny // 2) * dy
    X, Y = np.meshgrid(x, y, indexing="xy")
    return PeriodicGrid(
        nx=nx,
        ny=ny,
        lx=float((nx - 1) * dx),
        ly=float((ny - 1) * dy),
        x=x,
        y=y,
        X=X,
        Y=Y,
        center_x=0.0,
        center_y=0.0,
    )
