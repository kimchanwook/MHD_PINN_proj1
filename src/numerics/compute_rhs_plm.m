function rhs = compute_rhs_plm(U, gridData, gamma)
% compute_rhs_plm
%
% Computes the semi-discrete finite-volume right-hand side using
% piecewise-linear reconstruction with a minmod limiter and Rusanov fluxes.
%
% INPUTS:
%   U      - conserved state on physical domain, size [Ny, Nx, 8]
%   gridData   - gridData structure from make_uniform_grid
%   gamma  - ratio of specific heats
%
% OUTPUT:
%   rhs    - time derivative dU/dt on physical domain, size [Ny, Nx, 8]

Ny = gridData.Ny;
Nx = gridData.Nx;

[ULx, URx, ULy, URy] = reconstruct_plm_minmod(U);
Fx_num = rusanov_flux_x(ULx, URx, gamma);
Gy_num = rusanov_flux_y(ULy, URy, gamma);

rhs_x = -(Fx_num(:, 2:Nx+1, :) - Fx_num(:, 1:Nx,   :)) ./ gridData.dx;
rhs_y = -(Gy_num(2:Ny+1, :, :) - Gy_num(1:Ny,   :, :)) ./ gridData.dy;

rhs = rhs_x + rhs_y;

end
