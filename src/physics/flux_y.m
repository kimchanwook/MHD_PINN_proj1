function Gy = flux_y(U, gamma)
% flux_y
%
% Computes the physical flux vector in the y-direction for ideal MHD.

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

Gy = zeros(size(U));

Gy(:,:,1) = my;
Gy(:,:,2) = my .* ux - By .* Bx;
Gy(:,:,3) = my .* uy + ptot - By.^2;
Gy(:,:,4) = my .* uz - By .* Bz;
Gy(:,:,5) = (E + ptot) .* uy - udotB .* By;
Gy(:,:,6) = uy .* Bx - ux .* By;
Gy(:,:,7) = 0;
Gy(:,:,8) = uy .* Bz - uz .* By;

end
