"""Small finite-difference/FFT tools for periodic vorticity experiments."""

from .diagnostics import circulation, enstrophy, vorticity_extrema
from .dynamics import vorticity_tendency
from .evolution import EvolutionConfig, EvolutionSnapshot, evolve_state
from .grid import PeriodicGrid, make_centered_periodic_grid, make_periodic_grid
from .jacobian import arakawa_jacobian, centered_jacobian, jacobian
from .poisson import solve_poisson_fft
from .timestepping import rk4_step_state, rk4_step_vorticity

__all__ = [
    "PeriodicGrid",
    "EvolutionConfig",
    "EvolutionSnapshot",
    "arakawa_jacobian",
    "centered_jacobian",
    "circulation",
    "enstrophy",
    "evolve_state",
    "jacobian",
    "make_centered_periodic_grid",
    "make_periodic_grid",
    "rk4_step_state",
    "rk4_step_vorticity",
    "solve_poisson_fft",
    "vorticity_extrema",
    "vorticity_tendency",
]
