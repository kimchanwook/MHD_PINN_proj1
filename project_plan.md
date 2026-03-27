# Project plan

## Project title
Physics-Informed Neural Network for Two-Dimensional MHD Evolution

## Test cases
1. Alfvén wave propagation for verification
2. Magnetic pressure pulse for nonlinear demonstration

## Development order
1. Build conventional 2D ideal MHD solver
2. Verify on small-amplitude Alfvén wave
3. Demonstrate nonlinear behavior with magnetic pressure pulse
4. Add higher-order numerics and stronger diagnostics
5. Build PINN constrained by the ideal MHD equations
6. Compare PINN predictions against the conventional solver

## Current milestone status
### Completed so far
- Folder structure
- Primitive/conserved conversion routines
- Pressure and total-energy helper functions
- Physical flux functions
- Grid generation
- Fast-wave-speed routines
- CFL time-step routine
- Periodic ghost-cell boundary condition routine
- First-order Rusanov flux routines
- First-order right-hand-side assembly
- Forward Euler update
- Uniform-state preservation test
- Initial Alfvén-wave setup and case runner
- L2 error diagnostic
- Phase-speed estimate
- Line-cut plotting for Alfvén verification

### Next planned items
- RK2 time stepping
- Divergence diagnostic
- Convergence study
- Piecewise-linear reconstruction with limiter
- Magnetic pressure pulse case
- Physics notes as .tex + .pdf
- PINN formulation and comparison
