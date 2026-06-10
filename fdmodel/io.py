"""Read and write model state output files."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from .fields import ModelState
from .grid import make_periodic_grid


def save_state_npz(
    path: str | Path,
    state: ModelState,
    *,
    case: str,
    component: str,
    metadata: dict[str, Any] | None = None,
) -> None:
    """Save grid coordinates and model state arrays in a compressed NPZ file."""

    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    metadata = metadata or {}

    np.savez_compressed(
        output_path,
        case=np.array(case),
        component=np.array(component),
        nx=np.array(state.grid.nx),
        ny=np.array(state.grid.ny),
        lx=np.array(state.grid.lx),
        ly=np.array(state.grid.ly),
        x=state.grid.x,
        y=state.grid.y,
        zeta=state.zeta,
        psi=state.psi,
        u=state.u,
        v=state.v,
        metadata_keys=np.array(list(metadata.keys()), dtype=str),
        metadata_values=np.array([str(value) for value in metadata.values()], dtype=str),
    )


def load_state_npz(path: str | Path) -> ModelState:
    """Load a state file written by :func:`save_state_npz`."""

    data = np.load(Path(path), allow_pickle=False)
    grid = make_periodic_grid(
        int(data["nx"]),
        int(data["ny"]),
        lx=float(data["lx"]),
        ly=float(data["ly"]),
    )
    return ModelState(
        grid=grid,
        zeta=np.array(data["zeta"]),
        psi=np.array(data["psi"]),
        u=np.array(data["u"]),
        v=np.array(data["v"]),
    )


def load_metadata_npz(path: str | Path) -> dict[str, str]:
    """Load string metadata from a state file."""

    data = np.load(Path(path), allow_pickle=False)
    keys = [str(key) for key in data["metadata_keys"]]
    values = [str(value) for value in data["metadata_values"]]
    metadata = dict(zip(keys, values, strict=True))
    metadata["case"] = str(data["case"])
    metadata["component"] = str(data["component"])
    return metadata
