# Development Roadmap

## Phase 0

Project setup

- Initialize Git repository
- Create documentation
- Establish Python environment

Status: In Progress

---

## Phase 1

Static Model Infrastructure

Deliverables:

- grid generation
- monopole initialization
- perturbation initialization
- circulation correction
- FFT Poisson inversion
- wind diagnostics

Success Criteria:

- Initial vortex created
- Streamfunction successfully inverted
- Winds computed

---

## Phase 2

Core Dynamics

Deliverables:

- finite-difference operators
- Jacobian operator
- Laplacian diffusion
- RHS tendency calculation

Success Criteria:

- Stable tendency computation

---

## Phase 3

Time Integration

Deliverables:

- RK4 implementation
- model state updates

Success Criteria:

- Stable short integrations

---

## Phase 4

Monopole Experiment

Reproduce thesis axisymmetrization experiment.

Success Criteria:

- Qualitative agreement with original results

---

## Phase 5

Ring Vortex Experiment

Reproduce barotropic instability experiment.

Success Criteria:

- Polygonal eyewall evolution
- Vorticity mixing

---

## Phase 6

Validation

Compare against:

- thesis figures
- original model behavior
- conservation properties

---

## Phase 7

Modern Enhancements

Potential additions:

- Numba acceleration
- xarray support
- NetCDF output
- animation generation
- parameter sweeps