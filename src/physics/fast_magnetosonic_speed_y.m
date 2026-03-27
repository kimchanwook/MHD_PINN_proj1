function cfy = fast_magnetosonic_speed_y(U, gamma)
% fast_magnetosonic_speed_y
%
% Computes the fast magnetosonic speed in the y-direction for ideal MHD.

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
by2 = By.^2 ./ rho;

discriminant = (a2 + b2).^2 - 4 .* a2 .* by2;
discriminant = max(discriminant, 0);

cfy2 = 0.5 .* (a2 + b2 + sqrt(discriminant));
cfy  = sqrt(cfy2);

end
