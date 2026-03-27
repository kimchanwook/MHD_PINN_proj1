function cfx = fast_magnetosonic_speed_x(U, gamma)
% fast_magnetosonic_speed_x
%
% Computes the fast magnetosonic speed in the x-direction for ideal MHD.

rho = U(:,:,1);
mx  = U(:,:,2);
my  = U(:,:,3);
mz  = U(:,:,4);
E   = U(:,:,5);
Bx  = U(:,:,6);
By  = U(:,:,7);
Bz  = U(:,:,8);

p = compute_pressure(rho, mx, my, mz, E, Bx, By, Bz, gamma);

a2  = gamma .* p ./ rho;
b2  = (Bx.^2 + By.^2 + Bz.^2) ./ rho;
bx2 = Bx.^2 ./ rho;

discriminant = (a2 + b2).^2 - 4 .* a2 .* bx2;
discriminant = max(discriminant, 0);

cfx2 = 0.5 .* (a2 + b2 + sqrt(discriminant));
cfx  = sqrt(cfx2);

end
