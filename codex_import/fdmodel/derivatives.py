"""Periodic centered finite-difference derivatives."""

from __future__ import annotations

import numpy as np

from .grid import PeriodicGrid


def d_dx(field: np.ndarray, grid: PeriodicGrid) -> np.ndarray:
    """Centered periodic derivative with respect to x."""

    if field.shape != grid.shape:
        raise ValueError(f"field shape {field.shape} does not match grid shape {grid.shape}")
    return (np.roll(field, -1, axis=1) - np.roll(field, 1, axis=1)) / (2.0 * grid.dx)


def d_dy(field: np.ndarray, grid: PeriodicGrid) -> np.ndarray:
    """Centered periodic derivative with respect to y."""

    if field.shape != grid.shape:
        raise ValueError(f"field shape {field.shape} does not match grid shape {grid.shape}")
    return (np.roll(field, -1, axis=0) - np.roll(field, 1, axis=0)) / (2.0 * grid.dy)


def laplacian(field: np.ndarray, grid: PeriodicGrid) -> np.ndarray:
    """Second-order centered periodic Laplacian."""

    if field.shape != grid.shape:
        raise ValueError(f"field shape {field.shape} does not match grid shape {grid.shape}")
    d2x = (np.roll(field, -1, axis=1) - 2.0 * field + np.roll(field, 1, axis=1)) / (grid.dx * grid.dx)
    d2y = (np.roll(field, -1, axis=0) - 2.0 * field + np.roll(field, 1, axis=0)) / (grid.dy * grid.dy)
    return d2x + d2y
