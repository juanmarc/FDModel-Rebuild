"""Multi-step evolution control for vorticity simulations."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass

from .fields import ModelState
from .timestepping import rk4_step_state


@dataclass(frozen=True)
class EvolutionConfig:
    """Time-control parameters for a fixed-step integration."""

    dt: float
    viscosity: float
    max_time_count: int
    output_interval: int
    jacobian_method: str = "arakawa"

    def __post_init__(self) -> None:
        if self.dt <= 0.0:
            raise ValueError("dt must be positive")
        if self.viscosity < 0.0:
            raise ValueError("viscosity must be nonnegative")
        if self.max_time_count < 0:
            raise ValueError("max_time_count must be nonnegative")
        if self.output_interval <= 0:
            raise ValueError("output_interval must be positive")


@dataclass(frozen=True)
class EvolutionSnapshot:
    """State output at one model time."""

    step: int
    time_seconds: float
    state: ModelState


def evolve_state(initial_state: ModelState, config: EvolutionConfig) -> Iterator[EvolutionSnapshot]:
    """Yield model states at step 0 and each configured output interval."""

    state = initial_state
    yield EvolutionSnapshot(step=0, time_seconds=0.0, state=state)

    for step in range(1, config.max_time_count + 1):
        state = rk4_step_state(
            state,
            dt=config.dt,
            viscosity=config.viscosity,
            jacobian_method=config.jacobian_method,
        )
        if step % config.output_interval == 0 or step == config.max_time_count:
            yield EvolutionSnapshot(step=step, time_seconds=step * config.dt, state=state)
