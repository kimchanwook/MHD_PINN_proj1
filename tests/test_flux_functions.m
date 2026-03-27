function test_flux_functions()
% test_flux_functions
%
% Basic sanity tests for ideal MHD flux functions in x and y.

clc;

gamma = 5/3;

Ny = 3;
Nx = 4;
nVar = 8;

V = zeros(Ny, Nx, nVar);

rho0 = 1.0;
ux0  = 0.2;
uy0  = -0.1;
uz0  = 0.05;
p0   = 0.8;
Bx0  = 0.7;
By0  = -0.3;
Bz0  = 0.1;

V(:,:,1) = rho0;
V(:,:,2) = ux0;
V(:,:,3) = uy0;
V(:,:,4) = uz0;
V(:,:,5) = p0;
V(:,:,6) = Bx0;
V(:,:,7) = By0;
V(:,:,8) = Bz0;

U  = primitive_to_conserved(V, gamma);
Fx = flux_x(U, gamma);
Gy = flux_y(U, gamma);

assert(isequal(size(Fx), size(U)), 'FAIL: flux_x size mismatch.');
assert(isequal(size(Gy), size(U)), 'FAIL: flux_y size mismatch.');

mx = U(:,:,2);
my = U(:,:,3);

tmp = Fx(:,:,1) - mx;
err_mass_x = max(abs(tmp(:)));
tmp = Gy(:,:,1) - my;
err_mass_y = max(abs(tmp(:)));

assert(err_mass_x < 1e-12, 'FAIL: x mass flux is not equal to mx.');
assert(err_mass_y < 1e-12, 'FAIL: y mass flux is not equal to my.');

tmp = Fx(:,:,6);
err_Bx_flux_x = max(abs(tmp(:)));
assert(err_Bx_flux_x < 1e-12, 'FAIL: x-direction flux of Bx is not zero.');

tmp = Gy(:,:,7);
err_By_flux_y = max(abs(tmp(:)));
assert(err_By_flux_y < 1e-12, 'FAIL: y-direction flux of By is not zero.');

ptot0 = p0 + 0.5 * (Bx0^2 + By0^2 + Bz0^2);
expected_Fx2 = rho0 * ux0^2 + ptot0 - Bx0^2;

tmp = Fx(:,:,2) - expected_Fx2;
err_Fx2 = max(abs(tmp(:)));
assert(err_Fx2 < 1e-12, 'FAIL: x-momentum flux does not match expected value.');

expected_Gy3 = rho0 * uy0^2 + ptot0 - By0^2;
tmp = Gy(:,:,3) - expected_Gy3;
err_Gy3 = max(abs(tmp(:)));
assert(err_Gy3 < 1e-12, 'FAIL: y-momentum flux does not match expected value.');

fprintf('PASS: Flux function tests passed.\n');

end
