import numpy as np

from fdmodel.diagnostics import azimuthal_mean_profile, radial_wind_profile
from fdmodel.grid import make_periodic_grid
from fdmodel.initial_conditions import (
    monopole_initial_states,
    schubert_ring_initial_states,
)
from fdmodel.io import load_metadata_npz, load_state_npz, save_state_npz


def test_monopole_initial_states_include_base_perturbation_and_total():
    grid = make_periodic_grid(96, 96, lx=192.0e3, ly=192.0e3)

    states = monopole_initial_states(grid)

    assert set(states) == {"base", "perturbation", "total"}
    for state in states.values():
        assert state.zeta.shape == grid.shape
        assert state.psi.shape == grid.shape
        assert state.u.shape == grid.shape
        assert state.v.shape == grid.shape
        assert abs(np.mean(state.zeta)) < 1.0e-18
        assert abs(np.mean(state.psi)) < 1.0e-14


def test_ring_initial_states_include_base_perturbation_and_total():
    grid = make_periodic_grid(101, 101, lx=202.0e3, ly=202.0e3)

    states = schubert_ring_initial_states(grid)

    assert set(states) == {"base", "perturbation", "total"}
    assert states["base"].zeta.max() > states["base"].zeta[0, 0]
    assert states["perturbation"].zeta.max() > 0.0
    assert states["perturbation"].zeta.min() < 0.0


def test_state_npz_round_trip_preserves_arrays(tmp_path):
    grid = make_periodic_grid(48, 48, lx=96.0e3, ly=96.0e3)
    state = monopole_initial_states(grid)["base"]
    path = tmp_path / "monopole_base.npz"

    save_state_npz(path, state, case="monopole", component="base", metadata={"dx": grid.dx})
    loaded = load_state_npz(path)
    metadata = load_metadata_npz(path)

    assert loaded.grid.shape == state.grid.shape
    assert np.allclose(loaded.zeta, state.zeta)
    assert np.allclose(loaded.psi, state.psi)
    assert np.allclose(loaded.u, state.u)
    assert np.allclose(loaded.v, state.v)
    assert metadata["case"] == "monopole"
    assert metadata["component"] == "base"


def test_radial_profiles_return_expected_lengths():
    grid = make_periodic_grid(64, 64, lx=128.0e3, ly=128.0e3)
    state = monopole_initial_states(grid)["base"]

    radius, zeta_mean, counts = azimuthal_mean_profile(
        state.zeta,
        grid,
        bin_width=2.0e3,
        max_radius=40.0e3,
    )
    wind_radius, tangential_wind, wind_counts = radial_wind_profile(
        state.u,
        state.v,
        grid,
        bin_width=2.0e3,
        max_radius=40.0e3,
    )

    assert radius.shape == zeta_mean.shape == counts.shape
    assert wind_radius.shape == tangential_wind.shape == wind_counts.shape
    assert np.nanmax(zeta_mean) > np.nanmin(zeta_mean)
