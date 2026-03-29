function totalMass = compute_total_mass(U, gridData)
% compute_total_mass
%
% Computes the total mass in the computational domain.

rho = U(:,:,1);
totalMass = sum(rho(:)) * gridData.dx * gridData.dy;

end
