import numpy as np

from fdmodel.diagnostics import azimuthal_mean_profile, radial_wind_profile
from fdmodel.grid import make_centered_periodic_grid, make_periodic_grid
from fdmodel.initial_conditions import (
    monopole_initial_states,
    schubert_ring_initial_states,
)
from fdmodel.io import load_metadata_npz, load_state_npz, save_state_npz


def trapezoidal_circulation(zeta, grid):
    return np.trapezoid(np.trapezoid(zeta, x=grid.x, axis=1), x=grid.y, axis=0)


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


def test_monopole_base_state_uses_centered_grid_origin():
    grid = make_centered_periodic_grid(301, 301, dx=2.0e3, dy=2.0e3)

    states = monopole_initial_states(grid)
    max_index = np.unravel_index(np.argmax(states["base"].zeta), grid.shape)

    assert max_index == (150, 150)
    assert np.isclose(grid.x[max_index[1]], 0.0)
    assert np.isclose(grid.y[max_index[0]], 0.0)
    assert abs(trapezoidal_circulation(states["base"].zeta, grid)) < 1.0e-7


def test_ring_initial_states_include_base_perturbation_and_total():
    grid = make_periodic_grid(101, 101, lx=202.0e3, ly=202.0e3)

    states = schubert_ring_initial_states(grid)

    assert set(states) == {"base", "perturbation", "total"}
    assert states["base"].zeta.max() > states["base"].zeta[0, 0]
    assert states["perturbation"].zeta.max() > 0.0
    assert states["perturbation"].zeta.min() < 0.0


def test_ring_initial_states_use_centered_grid_origin():
    grid = make_centered_periodic_grid(301, 301, dx=2.0e3, dy=2.0e3)

    states = schubert_ring_initial_states(grid)

    assert states["base"].zeta.shape == grid.shape
    assert abs(trapezoidal_circulation(states["total"].zeta, grid)) < 1.0e-7
    assert states["base"].zeta[150, 150] > 0.0
    assert states["base"].zeta[0, 0] < 0.0


def test_ring_perturbation_uses_eight_percent_peak_amplitude():
    grid = make_centered_periodic_grid(301, 301, dx=2.0e3, dy=2.0e3)

    perturbation = schubert_ring_initial_states(grid)["perturbation"].zeta

    assert np.isclose(float(np.max(perturbation)), 5.6e-4, rtol=2.0e-2)
    assert float(np.min(perturbation)) < -4.0e-5


def test_state_npz_round_trip_preserves_arrays(tmp_path):
    grid = make_periodic_grid(48, 48, lx=96.0e3, ly=96.0e3)
    state = monopole_initial_states(grid)["base"]
    path = tmp_path / "monopole_base.npz"

    save_state_npz(path, state, case="monopole", component="base", metadata={"dx": grid.dx})
    loaded = load_state_npz(path)
    metadata = load_metadata_npz(path)

    assert loaded.grid.shape == state.grid.shape
    assert loaded.grid.center == state.grid.center
    assert np.allclose(loaded.grid.x, state.grid.x)
    assert np.allclose(loaded.grid.y, state.grid.y)
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
