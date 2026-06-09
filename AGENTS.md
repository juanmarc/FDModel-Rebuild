# FDModel-Rebuild

## Repository Root

/Users/marc/projects/FDModel-Rebuild

## Project Purpose

This project is a modern Python rebuild of a finite-difference barotropic nondivergent atmospheric model originally developed in Fortran by Marc Hidalgo as part of an M.S. thesis at Colorado State University.

The primary objective is to reproduce the original finite-difference model and associated hurricane vortex dynamics experiments using modern Python tools and software engineering practices.

## Instructions for Coding Agents

Before making any changes:

1. Read:
   - docs/PROJECT_CONTEXT.md
   - docs/ARCHITECTURE.md
   - docs/ROADMAP.md
   - docs/FORTRAN_INVENTORY.md

2. Preserve the existing project structure.

3. Do not create new top-level directories without approval.

4. Python source code belongs under:
   - fdmodel/

5. Run scripts belong under:
   - runs/

6. Tests belong under:
   - tests/

7. Documentation belongs under:
   - docs/

8. Favor readability and scientific correctness over optimization.

9. Maintain modular design.

10. Add tests for new numerical functionality whenever practical.

## Current Development Phase

Phase 1 complete:

- periodic Cartesian grid
- monopole initialization
- perturbation initialization
- circulation correction
- FFT Poisson inversion
- velocity diagnostics
- basic diagnostics
- unit tests

Next phase:

- Jacobian operator
- Laplacian operator
- RHS tendency computation
- RK4 time integration

## Coding Standards

- Python 3.12+
- NumPy
- SciPy
- pytest

Use type hints where practical.

Prefer small, focused functions with clear docstrings.