"""Plot model output fields and radial diagnostics."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .diagnostics import azimuthal_mean_profile, radial_wind_profile
from .fields import ModelState
from .grid import PeriodicGrid
from .io import load_state_npz


def plot_vorticity_plan_view(
    state: ModelState,
    path: str | Path,
    *,
    title: str,
    extent_radius: float | None = None,
    levels: np.ndarray | None = None,
) -> None:
    """Plot a plan view of vorticity and save it to ``path``."""

    plt = _import_pyplot()
    grid = state.grid
    coord_scale = _coordinate_scale(grid)
    coord_label = "km" if coord_scale == 1000.0 else "model units"
    x = (grid.x - 0.5 * grid.lx) / coord_scale
    y = (grid.y - 0.5 * grid.ly) / coord_scale
    zeta = state.zeta

    if extent_radius is not None:
        extent = extent_radius / coord_scale
        xmask = np.abs(x) <= extent
        ymask = np.abs(y) <= extent
        x = x[xmask]
        y = y[ymask]
        zeta = zeta[np.ix_(ymask, xmask)]

    if levels is None:
        max_abs = float(np.nanmax(np.abs(zeta)))
        levels = np.linspace(-max_abs, max_abs, 17)

    fig, ax = plt.subplots(figsize=(6.0, 5.8), constrained_layout=True)
    contour_fill = ax.contourf(x, y, zeta, levels=levels, cmap="RdBu_r", extend="both")
    contours = ax.contour(x, y, zeta, levels=levels, colors="black", linewidths=0.6)
    ax.clabel(contours, inline=True, fontsize=8)
    fig.colorbar(contour_fill, ax=ax, label="vorticity")
    ax.set_title(title)
    ax.set_xlabel(f"x ({coord_label})")
    ax.set_ylabel(f"y ({coord_label})")
    ax.set_aspect("equal")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_radial_profiles(
    state: ModelState,
    path: str | Path,
    *,
    title: str,
    bin_width: float,
    max_radius: float | None = None,
) -> None:
    """Plot azimuthal-mean vorticity and tangential wind profiles."""

    plt = _import_pyplot()
    radius_zeta, zeta_mean, _ = azimuthal_mean_profile(
        state.zeta,
        state.grid,
        bin_width=bin_width,
        max_radius=max_radius,
    )
    radius_wind, tangential_wind, _ = radial_wind_profile(
        state.u,
        state.v,
        state.grid,
        bin_width=bin_width,
        max_radius=max_radius,
    )

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8), constrained_layout=True)
    axes[0].plot(radius_zeta, zeta_mean, color="black")
    axes[0].set_title("Azimuthal-mean vorticity")
    axes[0].set_xlabel("radius")
    axes[0].set_ylabel("vorticity")
    axes[0].grid(True, alpha=0.25)

    axes[1].plot(radius_wind, tangential_wind, color="black")
    axes[1].set_title("Azimuthal-mean tangential wind")
    axes[1].set_xlabel("radius")
    axes[1].set_ylabel("tangential wind")
    axes[1].grid(True, alpha=0.25)
    fig.suptitle(title)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_state_file(
    state_path: str | Path,
    output_dir: str | Path,
    *,
    title: str,
    bin_width: float,
    extent_radius: float | None = None,
    max_radius: float | None = None,
) -> tuple[Path, Path]:
    """Read one saved state and write plan/radial plots."""

    state = load_state_npz(state_path)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    stem = Path(state_path).stem
    plan_path = output / f"{stem}_init_plan.png"
    radial_path = output / f"{stem}_init_radial.png"
    plot_vorticity_plan_view(state, plan_path, title=title, extent_radius=extent_radius)
    plot_radial_profiles(state, radial_path, title=title, bin_width=bin_width, max_radius=max_radius)
    return plan_path, radial_path


def _import_pyplot():
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "Plotting requires matplotlib. Install it in the project environment "
            "before running plot commands."
        ) from exc
    return plt


def _coordinate_scale(grid: PeriodicGrid) -> float:
    return 1000.0 if max(grid.lx, grid.ly) > 10_000.0 else 1.0
