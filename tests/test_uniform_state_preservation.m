function test_uniform_state_preservation()
% test_uniform_state_preservation
%
% Tests whether a spatially uniform state is preserved by the first-order
% finite-volume update with periodic boundary conditions.

clc;

gamma = 5/3;
CFL   = 0.4;

Nx = 10;
Ny = 12;

gridData = make_uniform_grid(0.0, 1.0, 0.0, 1.0, Nx, Ny);

nVar = 8;
V = zeros(Ny, Nx, nVar);

rho0 = 1.1;
ux0  = 0.2;
uy0  = -0.15;
uz0  = 0.05;
p0   = 0.9;
Bx0  = 0.4;
By0  = -0.25;
Bz0  = 0.1;

V(:,:,1) = rho0;
V(:,:,2) = ux0;
V(:,:,3) = uy0;
V(:,:,4) = uz0;
V(:,:,5) = p0;
V(:,:,6) = Bx0;
V(:,:,7) = By0;
V(:,:,8) = Bz0;

U = primitive_to_conserved(V, gamma);

rhs = compute_rhs(U, gridData, gamma);
rhs_max = max(abs(rhs(:)));
assert(rhs_max < 1e-12, 'FAIL: Uniform-state rhs is not zero.');

dt = compute_time_step(U, gridData, gamma, CFL);
Unew = update_fv_euler(U, gridData, gamma, dt);

update_err = max(abs(Unew(:) - U(:)));
assert(update_err < 1e-12, 'FAIL: Uniform state changed after Euler update.');

fprintf('PASS: Uniform state is preserved under periodic first-order update.\n');

end
