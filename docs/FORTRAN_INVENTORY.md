# FORTRAN Inventory

## Purpose

This document catalogs the original Fortran code base used for the finite-difference
barotropic nondivergent vortex model and maps legacy functionality to the new Python
implementation.

The objective is to:

1. Understand the original code structure.
2. Preserve scientific intent.
3. Identify numerical methods used.
4. Track progress during the Python port.
5. Document any deviations from the original implementation.

---

# Inventory Status

| File | Purpose Identified | Python Replacement | Port Status |
|--------|--------|--------|--------|
| TBD | TBD | TBD | Not Started |

---

# Original Model Architecture

## High-Level Workflow

Current understanding:

1. Read model configuration.
2. Initialize grid.
3. Initialize vorticity field.
4. Apply circulation correction.
5. Invert vorticity to streamfunction.
6. Compute winds.
7. Compute tendencies.
8. Advance solution in time.
9. Compute diagnostics.
10. Write output.

This workflow will be verified against the original code.

---

# File Inventory

---

## File: ______________________

### Purpose

Description:

TBD

### Key Routines

| Routine | Description |
|----------|----------|
| TBD | TBD |

### Inputs

- TBD

### Outputs

- TBD

### Numerical Methods

- TBD

### Python Mapping

| Legacy Function | Python Module |
|-----------------|---------------|
| TBD | TBD |

### Notes

TBD

---

## File: ______________________

### Purpose

Description:

TBD

### Key Routines

| Routine | Description |
|----------|----------|
| TBD | TBD |

### Inputs

- TBD

### Outputs

- TBD

### Numerical Methods

- TBD

### Python Mapping

| Legacy Function | Python Module |
|-----------------|---------------|
| TBD | TBD |

### Notes

TBD

---

# Numerical Components Inventory

## Grid Generation

Status: Not Yet Reviewed

Relevant Files:

- TBD

Python Target:

- grid.py

---

## Initial Conditions

Status: Not Yet Reviewed

Relevant Files:

- TBD

Python Target:

- initial_conditions.py

---

## Poisson Solver

Status: Not Yet Reviewed

Relevant Files:

- TBD

Python Target:

- poisson.py

Notes:

Original code likely uses NCAR elliptic solver routines.

Current Python implementation will use FFT-based inversion.

---

## Finite Difference Operators

Status: Not Yet Reviewed

Relevant Files:

- TBD

Python Target:

- derivatives.py

---

## Jacobian Calculation

Status: Not Yet Reviewed

Relevant Files:

- TBD

Python Target:

- jacobian.py

Questions:

- Was Arakawa used?
- Was centered differencing used?
- Were alternative schemes available?

---

## Time Integration

Status: Not Yet Reviewed

Relevant Files:

- TBD

Python Target:

- timestepping.py

Expected Method:

- RK4

Reference:

Thesis Chapter 2.

---

## Diagnostics

Status: Not Yet Reviewed

Relevant Files:

- TBD

Python Target:

- diagnostics.py

Expected Diagnostics:

- circulation
- kinetic energy
- angular momentum
- enstrophy
- palinstrophy

---

## Output

Status: Not Yet Reviewed

Relevant Files:

- TBD

Python Target:

- io.py

---

# Porting Decisions Log

## Decision 001

Original:

Unknown

Replacement:

FFT Poisson inversion using NumPy/SciPy.

Reason:

Modern hardware and periodic domain make FFT inversion simpler and faster than
porting legacy NCAR elliptic solver code.

Status:

Approved

---

## Decision 002

Original:

Fortran 77/90 implementation.

Replacement:

Python 3.12+ with NumPy/SciPy.

Reason:

Maintainability and developer productivity.

Status:

Approved

---

# Open Questions

1. Which files are specific to the FD model?
2. Which files belong only to the semi-spectral model?
3. Which files are support utilities?
4. Are there any known bugs in the original code?
5. Which original experiments should be used for validation?