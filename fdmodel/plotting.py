"""Plot model output fields and radial diagnostics."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .diagnostics import azimuthal_mean_profile, radial_wind_profile
from .fields import ModelState
from .grid import PeriodicGrid
from .io import load_metadata_npz, load_state_npz


@dataclass(frozen=True)
class PlotConvention:
    contour_levels: np.ndarray
    gray_min: float | None
    gray_max: float | None
    black_min: float | None
    negative_contour_levels: np.ndarray | None = None


def plot_vorticity_plan_view(
    state: ModelState,
    path: str | Path,
    *,
    title: str,
    extent_radius: float | None = None,
    levels: np.ndarray | None = None,
    convention: PlotConvention | None = None,
    annotate_extrema: bool = True,
) -> None:
    """Plot a plan view of vorticity and save it to ``path``."""

    plt = _import_pyplot()
    grid = state.grid
    coord_scale = _coordinate_scale(grid)
    coord_label = "km" if coord_scale == 1000.0 else "model units"
    x = (grid.x - grid.center[0]) / coord_scale
    y = (grid.y - grid.center[1]) / coord_scale
    zeta = state.zeta

    if extent_radius is not None:
        extent = extent_radius / coord_scale
        xmask = np.abs(x) <= extent
        ymask = np.abs(y) <= extent
        x = x[xmask]
        y = y[ymask]
        zeta = zeta[np.ix_(ymask, xmask)]

    if convention is None:
        convention = _monopole_plot_convention()
    if levels is None:
        levels = convention.contour_levels

    fig, ax = plt.subplots(figsize=(6.0, 5.8), constrained_layout=True)
    ax.set_facecolor("white")
    _shade_vorticity_bands(ax, x, y, zeta, convention)
    negative_levels = convention.negative_contour_levels
    if negative_levels is not None:
        negative_levels = negative_levels[
            (negative_levels >= float(np.nanmin(zeta))) & (negative_levels <= float(np.nanmax(zeta)))
        ]
        negative_levels = np.sort(negative_levels)
        if negative_levels.size:
            negative_contours = ax.contour(
                x,
                y,
                zeta,
                levels=negative_levels,
                colors="black",
                linewidths=0.7,
                linestyles="dotted",
            )
            ax.clabel(negative_contours, inline=True, fontsize=8, fmt=_format_scientific_label)
    contour_levels = levels[(levels >= float(np.nanmin(zeta))) & (levels <= float(np.nanmax(zeta)))]
    if contour_levels.size:
        contours = ax.contour(
            x,
            y,
            zeta,
            levels=contour_levels,
            colors="black",
            linewidths=0.7,
            linestyles="solid",
        )
        ax.clabel(contours, inline=True, fontsize=8, fmt=_format_scientific_label)
    ax.set_title(title)
    ax.set_xlabel(f"x ({coord_label})")
    ax.set_ylabel(f"y ({coord_label})")
    ax.set_aspect("equal")
    if annotate_extrema:
        _annotate_extrema(ax, zeta)
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
    title: str | None = None,
    bin_width: float,
    extent_radius: float | None = None,
    max_radius: float | None = None,
) -> tuple[Path, Path]:
    """Read one saved state and write plan/radial plots."""

    state = load_state_npz(state_path)
    metadata = load_metadata_npz(state_path)
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    stem = Path(state_path).stem
    suffix = "init" if _is_initial_condition(metadata) else "evolution"
    plan_path = output / f"{stem}_{suffix}_plan.png"
    radial_path = output / f"{stem}_{suffix}_radial.png"
    if title is None:
        title = _title_from_metadata(metadata)
    convention = _plot_convention_from_metadata(metadata)
    plot_vorticity_plan_view(state, plan_path, title=title, extent_radius=extent_radius, convention=convention)
    plot_radial_profiles(state, radial_path, title=title, bin_width=bin_width, max_radius=max_radius)
    return plan_path, radial_path


def _title_from_metadata(metadata: dict[str, str]) -> str:
    if "time_hours" in metadata:
        return f"Vorticity at t={float(metadata['time_hours']):.2f} h"
    if metadata.get("component") == "evolution":
        return "Vorticity at t=0.00 h"
    case = metadata.get("case", "state")
    component = metadata.get("component", "")
    return f"{case} {component} initial conditions".strip()


def _is_initial_condition(metadata: dict[str, str]) -> bool:
    return metadata.get("component") != "evolution"


def _annotate_extrema(ax, zeta: np.ndarray) -> None:
    max_text = f"Max={float(np.nanmax(zeta)):.5g} s$^{{-1}}$"
    min_text = f"Min={float(np.nanmin(zeta)):.5g} s$^{{-1}}$"
    text_style = {
        "transform": ax.transAxes,
        "fontsize": 7.5,
        "verticalalignment": "bottom",
        "bbox": {"boxstyle": "square,pad=0.2", "facecolor": "white", "edgecolor": "none", "alpha": 0.75},
    }
    ax.text(0.02, 0.03, max_text, horizontalalignment="left", **text_style)
    ax.text(0.98, 0.03, min_text, horizontalalignment="right", **text_style)


def _format_scientific_label(value: float) -> str:
    coefficient, exponent = f"{value:.2e}".split("e")
    return f"{coefficient}x10$^{{{int(exponent)}}}$"


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


def _plot_convention_from_metadata(metadata: dict[str, str]) -> PlotConvention:
    if metadata.get("case") == "ring":
        if metadata.get("component") == "perturbation":
            return _ring_perturbation_plot_convention()
        return _ring_plot_convention()
    return _monopole_plot_convention()


def _monopole_plot_convention() -> PlotConvention:
    return PlotConvention(
        contour_levels=np.array([5.0e-5] + [level * 1.0e-4 for level in range(1, 11)], dtype=float),
        gray_min=5.0e-5,
        gray_max=1.0e-4,
        black_min=1.0e-3,
    )


def _ring_plot_convention() -> PlotConvention:
    return PlotConvention(
        contour_levels=np.array([2.0e-4] + [level * 1.0e-3 for level in range(1, 11)], dtype=float),
        gray_min=2.0e-4,
        gray_max=1.0e-3,
        black_min=6.0e-3,
    )


def _ring_perturbation_plot_convention() -> PlotConvention:
    return PlotConvention(
        contour_levels=np.arange(7.5e-5, 15.0e-4 + 7.5e-5, 7.5e-5, dtype=float),
        gray_min=None,
        gray_max=None,
        black_min=None,
        negative_contour_levels=-np.arange(4.0e-5, 8.0e-4 + 4.0e-5, 4.0e-5, dtype=float),
    )


def _shade_vorticity_bands(
    ax,
    x: np.ndarray,
    y: np.ndarray,
    zeta: np.ndarray,
    convention: PlotConvention,
) -> None:
    if convention.gray_min is None or convention.gray_max is None or convention.black_min is None:
        return
    zeta_min = float(np.nanmin(zeta))
    zeta_max = float(np.nanmax(zeta))
    if zeta_min < convention.gray_min and zeta_max > convention.gray_min:
        upper = min(convention.gray_max, zeta_max)
        if upper > convention.gray_min:
            ax.contourf(x, y, zeta, levels=[convention.gray_min, upper], colors=["0.6"])
    if zeta_max >= convention.black_min:
        ax.contourf(x, y, zeta, levels=[convention.black_min, zeta_max + 1.0e-12], colors=["black"])
