function totalEnergy = compute_total_energy_domain(U, grid)
% compute_total_energy_domain
%
% Computes the total energy integrated over the computational domain.

E = U(:,:,5);
totalEnergy = sum(E(:)) * grid.dx * grid.dy;

end
