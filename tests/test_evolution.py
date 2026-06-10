import pytest

from fdmodel.evolution import EvolutionConfig, evolve_state
from fdmodel.fields import state_from_vorticity
from fdmodel.grid import make_periodic_grid


def test_evolve_state_yields_initial_and_output_interval_snapshots():
    grid = make_periodic_grid(24, 24)
    state = state_from_vorticity(0.01 * (grid.X * 0.0), grid)
    config = EvolutionConfig(dt=1.0, viscosity=0.0, max_time_count=5, output_interval=2)

    snapshots = list(evolve_state(state, config))

    assert [snapshot.step for snapshot in snapshots] == [0, 2, 4, 5]
    assert [snapshot.time_seconds for snapshot in snapshots] == [0.0, 2.0, 4.0, 5.0]


def test_evolution_config_rejects_invalid_values():
    with pytest.raises(ValueError, match="dt"):
        EvolutionConfig(dt=0.0, viscosity=0.0, max_time_count=1, output_interval=1)
    with pytest.raises(ValueError, match="viscosity"):
        EvolutionConfig(dt=1.0, viscosity=-1.0, max_time_count=1, output_interval=1)
    with pytest.raises(ValueError, match="max_time_count"):
        EvolutionConfig(dt=1.0, viscosity=0.0, max_time_count=-1, output_interval=1)
    with pytest.raises(ValueError, match="output_interval"):
        EvolutionConfig(dt=1.0, viscosity=0.0, max_time_count=1, output_interval=0)
