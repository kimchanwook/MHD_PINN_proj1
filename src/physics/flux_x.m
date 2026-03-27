function Fx = flux_x(U, gamma)
% flux_x
%
% Computes the physical flux vector in the x-direction for ideal MHD.

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

Bsq   = Bx.^2 + By.^2 + Bz.^2;
ptot  = p + 0.5 .* Bsq;
udotB = ux .* Bx + uy .* By + uz .* Bz;

Fx = zeros(size(U));

Fx(:,:,1) = mx;
Fx(:,:,2) = mx .* ux + ptot - Bx.^2;
Fx(:,:,3) = mx .* uy - Bx .* By;
Fx(:,:,4) = mx .* uz - Bx .* Bz;
Fx(:,:,5) = (E + ptot) .* ux - udotB .* Bx;
Fx(:,:,6) = 0;
Fx(:,:,7) = ux .* By - uy .* Bx;
Fx(:,:,8) = ux .* Bz - uz .* Bx;

end
