function test_rk2_uniform_state_preservation()
% test_rk2_uniform_state_preservation
%
% Verifies that a spatially uniform state remains unchanged under one RK2
% update with periodic boundary conditions.
%
% PHYSICAL EXPECTATION:
%   If every cell contains exactly the same state, all interface jumps are
%   zero. Therefore, every numerical flux is uniform and its discrete
%   divergence vanishes. Hence dU/dt = 0 and the RK2 update must leave the
%   state unchanged up to roundoff.

clc;

gamma = 5/3;
CFL   = 0.4;

Nx = 14;
Ny = 11;

gridData = make_uniform_grid(0.0, 1.0, 0.0, 1.0, Nx, Ny);

nVar = 8;
V = zeros(Ny, Nx, nVar);

rho0 = 1.2;
ux0  = -0.18;
uy0  = 0.11;
uz0  = 0.04;
p0   = 0.95;
Bx0  = 0.55;
By0  = -0.20;
Bz0  = 0.07;

V(:,:,1) = rho0;
V(:,:,2) = ux0;
V(:,:,3) = uy0;
V(:,:,4) = uz0;
V(:,:,5) = p0;
V(:,:,6) = Bx0;
V(:,:,7) = By0;
V(:,:,8) = Bz0;

U = primitive_to_conserved(V, gamma);
dt = compute_time_step(U, gridData, gamma, CFL);

Unew = update_fv_rk2(U, gridData, gamma, dt);
updateErr = max(abs(Unew(:) - U(:)));

assert(updateErr < 1e-12, 'FAIL: Uniform state changed after RK2 update.');

fprintf('PASS: Uniform state is preserved under periodic RK2 update.\n');

end
