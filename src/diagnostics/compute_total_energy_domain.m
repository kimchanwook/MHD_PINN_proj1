function totalEnergy = compute_total_energy_domain(U, gridData)
% compute_total_energy_domain
%
% Computes the total energy integrated over the computational domain.

E = U(:,:,5);
totalEnergy = sum(E(:)) * gridData.dx * gridData.dy;

end
