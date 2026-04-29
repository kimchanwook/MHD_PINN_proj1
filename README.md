# Physics-Informed Neural Network for Two-Dimensional MHD Evolution

This repository contains the early-stage MATLAB implementation of a two-dimensional ideal magnetohydrodynamics (MHD) solver that will serve as the conventional solver backbone before adding a physics-informed neural network (PINN). `This work was conducted with the assistance of a large language model (LLM), specifically ChatGPT`.

## Current contents
- Primitive/conserved variable conversion
- Total energy and pressure recovery
- Ideal MHD physical fluxes in x and y
- Uniform Cartesian grid generation
- Fast magnetosonic speed estimates in x and y
- CFL-based time step selection
- First-order finite-volume update with Rusanov fluxes
- Periodic ghost-cell boundary conditions
- First Alfvén-wave initialization and verification case
- Basic diagnostics and plotting
- Unit tests for the above components

## Current limitations
- First-order spatial reconstruction only
- Forward Euler time stepping only
- No divergence-cleaning yet
- No RK2 yet
- No magnetic pressure pulse case yet
- No PINN yet

## Quick start
From the project root in MATLAB:

```matlab
addpath(genpath(pwd));

test_variable_conversion
test_flux_functions
test_time_step_and_wave_speeds
test_uniform_state_preservation

results = alfven_wave_case();
```


## Stage 2 additions

This stage adds:
- piecewise-linear reconstruction with a minmod limiter,
- second-order Runge-Kutta updates using the PLM spatial operator,
- a nonlinear magnetic-pressure pulse case,
- a second physics note covering the Alfv'en verification case and the pulse demonstration.
