function V = conserved_to_primitive(U, gamma)
% conserved_to_primitive
%
% Converts conserved variables to primitive variables for ideal MHD.
%
% INPUT:
%   U(:,:,1) = rho  = mass density
%   U(:,:,2) = mx   = momentum density in x-direction
%   U(:,:,3) = my   = momentum density in y-direction
%   U(:,:,4) = mz   = momentum density in z-direction
%   U(:,:,5) = E    = total energy density
%   U(:,:,6) = Bx   = magnetic field in x-direction
%   U(:,:,7) = By   = magnetic field in y-direction
%   U(:,:,8) = Bz   = magnetic field in z-direction
%
%   gamma    = ratio of specific heats
%
% OUTPUT:
%   V(:,:,1) = rho
%   V(:,:,2) = ux
%   V(:,:,3) = uy
%   V(:,:,4) = uz
%   V(:,:,5) = p
%   V(:,:,6) = Bx
%   V(:,:,7) = By
%   V(:,:,8) = Bz

rho = U(:,:,1);
mx  = U(:,:,2);
my  = U(:,:,3);
mz  = U(:,:,4);
E   = U(:,:,5);
Bx  = U(:,:,6);
By  = U(:,:,7);
Bz  = U(:,:,8);

ux = mx ./ rho;
uy = my ./ rho;
uz = mz ./ rho;

p = compute_pressure(rho, mx, my, mz, E, Bx, By, Bz, gamma);

V = zeros(size(U));

V(:,:,1) = rho;
V(:,:,2) = ux;
V(:,:,3) = uy;
V(:,:,4) = uz;
V(:,:,5) = p;
V(:,:,6) = Bx;
V(:,:,7) = By;
V(:,:,8) = Bz;

end
