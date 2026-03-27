function test_time_step_and_wave_speeds()
% test_time_step_and_wave_speeds
%
% Tests grid generation, fast magnetosonic speeds, and CFL time step.

clc;

gamma = 5/3;
CFL   = 0.4;

xMin = 0.0;
xMax = 1.0;
yMin = -0.5;
yMax = 0.5;
Nx   = 8;
Ny   = 6;

grid = make_uniform_grid(xMin, xMax, yMin, yMax, Nx, Ny);

assert(abs(grid.dx - (xMax - xMin)/Nx) < 1e-14, 'FAIL: Incorrect dx.');
assert(abs(grid.dy - (yMax - yMin)/Ny) < 1e-14, 'FAIL: Incorrect dy.');
assert(isequal(size(grid.Xc), [Ny, Nx]), 'FAIL: Incorrect Xc size.');
assert(isequal(size(grid.Yc), [Ny, Nx]), 'FAIL: Incorrect Yc size.');

nVar = 8;
V = zeros(Ny, Nx, nVar);

rho0 = 1.0;
ux0  = 0.3;
uy0  = -0.2;
uz0  = 0.0;
p0   = 1.2;
Bx0  = 0.5;
By0  = 0.1;
Bz0  = -0.2;

V(:,:,1) = rho0;
V(:,:,2) = ux0;
V(:,:,3) = uy0;
V(:,:,4) = uz0;
V(:,:,5) = p0;
V(:,:,6) = Bx0;
V(:,:,7) = By0;
V(:,:,8) = Bz0;

U = primitive_to_conserved(V, gamma);

cfx = fast_magnetosonic_speed_x(U, gamma);
cfy = fast_magnetosonic_speed_y(U, gamma);

assert(all(cfx(:) > 0), 'FAIL: cfx is not strictly positive.');
assert(all(cfy(:) > 0), 'FAIL: cfy is not strictly positive.');

dt = compute_time_step(U, grid, gamma, CFL);
assert(isfinite(dt), 'FAIL: dt is not finite.');
assert(dt > 0, 'FAIL: dt is not positive.');

V_hydro = zeros(Ny, Nx, nVar);

rho_h = 1.4;
ux_h  = 0.0;
uy_h  = 0.0;
uz_h  = 0.0;
p_h   = 0.9;
Bx_h  = 0.0;
By_h  = 0.0;
Bz_h  = 0.0;

V_hydro(:,:,1) = rho_h;
V_hydro(:,:,2) = ux_h;
V_hydro(:,:,3) = uy_h;
V_hydro(:,:,4) = uz_h;
V_hydro(:,:,5) = p_h;
V_hydro(:,:,6) = Bx_h;
V_hydro(:,:,7) = By_h;
V_hydro(:,:,8) = Bz_h;

U_hydro = primitive_to_conserved(V_hydro, gamma);

cfx_h = fast_magnetosonic_speed_x(U_hydro, gamma);
cfy_h = fast_magnetosonic_speed_y(U_hydro, gamma);

sound_speed = sqrt(gamma * p_h / rho_h);

err_cfx = max(abs(cfx_h(:) - sound_speed));
err_cfy = max(abs(cfy_h(:) - sound_speed));

assert(err_cfx < 1e-12, 'FAIL: cfx does not reduce to sound speed when B = 0.');
assert(err_cfy < 1e-12, 'FAIL: cfy does not reduce to sound speed when B = 0.');

fprintf('PASS: Grid, wave-speed, and CFL time-step tests passed.\n');

end
