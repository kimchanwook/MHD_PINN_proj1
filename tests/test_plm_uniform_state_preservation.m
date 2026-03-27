function test_plm_uniform_state_preservation()
% test_plm_uniform_state_preservation
%
% Tests whether a spatially uniform state is preserved by the PLM + minmod
% + RK2 update with periodic boundary conditions.

clc;

gamma = 5/3;
CFL   = 0.4;
Nx = 14;
Ny = 10;
grid = make_uniform_grid(0.0, 1.0, 0.0, 1.0, Nx, Ny);

nVar = 8;
V = zeros(Ny, Nx, nVar);
V(:,:,1) = 1.05;
V(:,:,2) = 0.2;
V(:,:,3) = -0.12;
V(:,:,4) = 0.03;
V(:,:,5) = 0.95;
V(:,:,6) = 0.5;
V(:,:,7) = -0.2;
V(:,:,8) = 0.15;

U = primitive_to_conserved(V, gamma);
rhs = compute_rhs_plm(U, grid, gamma);
assert(max(abs(rhs(:))) < 1e-12, 'FAIL: Uniform-state rhs is not zero for PLM.');

dt = compute_time_step(U, grid, gamma, CFL);
Unew = update_fv_rk2_plm(U, grid, gamma, dt);
assert(max(abs(Unew(:) - U(:))) < 1e-12, 'FAIL: Uniform state changed after PLM RK2 update.');

fprintf('PASS: Uniform state is preserved under PLM + RK2 update.\n');

end
