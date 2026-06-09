# FD Model Rebuild

## Overview

This project is a modern Python rebuild of a finite-difference barotropic nondivergent atmospheric model originally developed by Marc Hidalgo in 1998-1999 as part of a Master of Science thesis at Colorado State University.

The original model was written in Fortran and used to investigate hurricane vortex dynamics, including:

- Axisymmetrization of finite-amplitude vorticity perturbations
- Vortex Rossby wave dynamics
- Barotropic instability
- Polygonal eyewall formation
- Vorticity mixing processes in hurricane cores

The original research is documented in:

Hidalgo, J.M. (1999)
"A Semi-Spectral Numerical Method for Modeling the Vorticity Dynamics of the Near-Core Region of Hurricane-Like Vortices"

Although the thesis focused heavily on a semi-spectral cylindrical model, this rebuild effort will focus exclusively on the fully finite-difference Cartesian model.

---

## Project Goals

### Scientific Goals

1. Recreate the original finite-difference model in Python.
2. Reproduce the monopole axisymmetrization experiments.
3. Reproduce the ring-vortex instability experiments.
4. Reproduce the polygonal eyewall evolution described by Schubert et al. (1999).
5. Validate results against the original thesis.

### Technical Goals

1. Refresh hands-on software development skills.
2. Refresh numerical modeling skills.
3. Demonstrate continued technical currency.
4. Apply modern Python scientific-computing practices.
5. Maintain a modular architecture suitable for future extension.

---

## Governing Equations

Barotropic nondivergent vorticity equation:

dζ/dt + J(ψ,ζ) = ν∇²ζ

where:

ζ = relative vorticity

ψ = streamfunction

J(ψ,ζ) = Jacobian

ν = viscosity coefficient

Streamfunction inversion:

∇²ψ = ζ

Velocity components:

u = -∂ψ/∂y

v =  ∂ψ/∂x

---

## Numerical Assumptions

- Cartesian x-y grid
- Doubly periodic boundaries
- Finite-difference spatial derivatives
- Runge-Kutta time integration
- FFT-based Poisson inversion
- NumPy-based implementation
- Python 3.12+

---

## Long-Term Vision

The final product should reproduce the original thesis experiments while leveraging modern hardware and software practices.