function rhs = compute_rhs(U, gridData, gamma)
% compute_rhs
%
% Computes the semi-discrete finite-volume right-hand side for ideal MHD
% using first-order piecewise-constant states and Rusanov numerical fluxes.

Ny = gridData.Ny;
Nx = gridData.Nx;

Ubc = apply_periodic_bc(U);

ULx = Ubc(2:Ny+1, 1:Nx+1,   :);
URx = Ubc(2:Ny+1, 2:Nx+2,   :);

Fx_num = rusanov_flux_x(ULx, URx, gamma);

ULy = Ubc(1:Ny+1,   2:Nx+1, :);
URy = Ubc(2:Ny+2,   2:Nx+1, :);

Gy_num = rusanov_flux_y(ULy, URy, gamma);

rhs_x = -(Fx_num(:, 2:Nx+1, :) - Fx_num(:, 1:Nx,   :)) ./ gridData.dx;
rhs_y = -(Gy_num(2:Ny+1, :, :) - Gy_num(1:Ny,   :, :)) ./ gridData.dy;

rhs = rhs_x + rhs_y;

end
