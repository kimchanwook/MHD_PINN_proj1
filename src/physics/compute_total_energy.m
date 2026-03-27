function E = compute_total_energy(rho, ux, uy, uz, p, Bx, By, Bz, gamma)
% compute_total_energy
%
% Computes the total energy density E for ideal magnetohydrodynamics (MHD).
%
% INPUTS:
%   rho   - mass density
%   ux    - velocity component in x-direction
%   uy    - velocity component in y-direction
%   uz    - velocity component in z-direction
%   p     - gas pressure
%   Bx    - magnetic field component in x-direction
%   By    - magnetic field component in y-direction
%   Bz    - magnetic field component in z-direction
%   gamma - ratio of specific heats
%
% OUTPUT:
%   E     - total energy density
%
% FORMULA:
%   E = p/(gamma - 1)
%       + 0.5*rho*(ux^2 + uy^2 + uz^2)
%       + 0.5*(Bx^2 + By^2 + Bz^2)

internal_energy = p ./ (gamma - 1);
kinetic_energy  = 0.5 .* rho .* (ux.^2 + uy.^2 + uz.^2);
magnetic_energy = 0.5 .* (Bx.^2 + By.^2 + Bz.^2);

E = internal_energy + kinetic_energy + magnetic_energy;

end
