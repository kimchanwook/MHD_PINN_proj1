function p = compute_pressure(rho, mx, my, mz, E, Bx, By, Bz, gamma)
% compute_pressure
%
% Computes the gas pressure p from conserved variables for ideal MHD.
%
% INPUTS:
%   rho   - mass density
%   mx    - momentum density in x-direction = rho * ux
%   my    - momentum density in y-direction = rho * uy
%   mz    - momentum density in z-direction = rho * uz
%   E     - total energy density
%   Bx    - magnetic field component in x-direction
%   By    - magnetic field component in y-direction
%   Bz    - magnetic field component in z-direction
%   gamma - ratio of specific heats
%
% OUTPUT:
%   p     - gas pressure

ux = mx ./ rho;
uy = my ./ rho;
uz = mz ./ rho;

kinetic_energy  = 0.5 .* rho .* (ux.^2 + uy.^2 + uz.^2);
magnetic_energy = 0.5 .* (Bx.^2 + By.^2 + Bz.^2);

p = (gamma - 1) .* (E - kinetic_energy - magnetic_energy);

end
