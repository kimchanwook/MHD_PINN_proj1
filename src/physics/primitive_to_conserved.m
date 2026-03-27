function U = primitive_to_conserved(V, gamma)
% primitive_to_conserved
%
% Converts primitive variables to conserved variables for ideal MHD.
%
% INPUT:
%   V(:,:,1) = rho  = mass density
%   V(:,:,2) = ux   = velocity in x-direction
%   V(:,:,3) = uy   = velocity in y-direction
%   V(:,:,4) = uz   = velocity in z-direction
%   V(:,:,5) = p    = gas pressure
%   V(:,:,6) = Bx   = magnetic field in x-direction
%   V(:,:,7) = By   = magnetic field in y-direction
%   V(:,:,8) = Bz   = magnetic field in z-direction
%
%   gamma    = ratio of specific heats
%
% OUTPUT:
%   U(:,:,1) = rho
%   U(:,:,2) = rho*ux
%   U(:,:,3) = rho*uy
%   U(:,:,4) = rho*uz
%   U(:,:,5) = E
%   U(:,:,6) = Bx
%   U(:,:,7) = By
%   U(:,:,8) = Bz

rho = V(:,:,1);
ux  = V(:,:,2);
uy  = V(:,:,3);
uz  = V(:,:,4);
p   = V(:,:,5);
Bx  = V(:,:,6);
By  = V(:,:,7);
Bz  = V(:,:,8);

mx = rho .* ux;
my = rho .* uy;
mz = rho .* uz;

E = compute_total_energy(rho, ux, uy, uz, p, Bx, By, Bz, gamma);

U = zeros(size(V));

U(:,:,1) = rho;
U(:,:,2) = mx;
U(:,:,3) = my;
U(:,:,4) = mz;
U(:,:,5) = E;
U(:,:,6) = Bx;
U(:,:,7) = By;
U(:,:,8) = Bz;

end
