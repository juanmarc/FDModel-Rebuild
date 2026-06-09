# Architecture

## Design Philosophy

The architecture should loosely mirror the original Fortran module structure while taking advantage of modern Python packaging practices.

Each module should have a single responsibility.

---

## Package Structure

fdmodel/

    config.py

    grid.py

    fields.py

    initial_conditions.py

    poisson.py

    derivatives.py

    jacobian.py

    dynamics.py

    timestepping.py

    diagnostics.py

    io.py

runs/

    run_monopole.py

    run_ring_vortex.py

tests/

    test_poisson.py

    test_derivatives.py

    test_diagnostics.py

---

## Module Responsibilities

### config.py

Simulation parameters

- grid dimensions
- timestep
- viscosity
- runtime

### grid.py

Grid generation

- x coordinates
- y coordinates
- mesh creation
- periodic indexing helpers

### fields.py

Container for model state

- vorticity
- streamfunction
- velocity

### initial_conditions.py

Initial vortex definitions

- monopole
- Gaussian perturbation
- ring vortex

### poisson.py

Poisson inversion

∇²ψ = ζ

Implemented using FFT methods.

### derivatives.py

Finite-difference operators

- d/dx
- d/dy
- Laplacian

### jacobian.py

Advection operator

J(ψ,ζ)

Initial implementation may use centered differences.

Future implementation should include Arakawa Jacobian.

### dynamics.py

Computes RHS tendency

dζ/dt

### timestepping.py

RK4 integration

### diagnostics.py

Scientific diagnostics

- circulation
- energy
- enstrophy
- palinstrophy

### io.py

Output support

- NumPy
- NetCDF (future)

---

## Coding Standards

- Type hints preferred
- Functions should remain short
- Unit tests for core numerical operators
- Clear docstrings