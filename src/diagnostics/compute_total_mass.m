function totalMass = compute_total_mass(U, grid)
% compute_total_mass
%
% Computes the total mass in the computational domain.

rho = U(:,:,1);
totalMass = sum(rho(:)) * grid.dx * grid.dy;

end
