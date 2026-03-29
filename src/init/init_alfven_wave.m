function [U0, exact, params] = init_alfven_wave(gridData, gamma, userParams)
% init_alfven_wave
%
% Initializes a small-amplitude shear Alfvén wave on a 2D gridData.

rho0 = userParams.rho0;
p0   = userParams.p0;
B0   = userParams.B0;
A    = userParams.A;
mode = userParams.mode;

Ny = gridData.Ny;
Nx = gridData.Nx;
nVar = 8;

Lx = gridData.xMax - gridData.xMin;
k  = 2 * pi * mode / Lx;
vA = B0 / sqrt(rho0);
lambda = Lx / mode;
period = lambda / vA;

X = gridData.Xc;

V0 = zeros(Ny, Nx, nVar);

sinPhase = sin(k .* X);

rho = rho0 .* ones(Ny, Nx);
ux  = zeros(Ny, Nx);
uy  = zeros(Ny, Nx);
uz  = A .* sinPhase;
p   = p0 .* ones(Ny, Nx);
Bx  = B0 .* ones(Ny, Nx);
By  = zeros(Ny, Nx);
Bz  = -sqrt(rho0) .* A .* sinPhase;

V0(:,:,1) = rho;
V0(:,:,2) = ux;
V0(:,:,3) = uy;
V0(:,:,4) = uz;
V0(:,:,5) = p;
V0(:,:,6) = Bx;
V0(:,:,7) = By;
V0(:,:,8) = Bz;

U0 = primitive_to_conserved(V0, gamma);

exact = struct();
exact.rho0   = rho0;
exact.p0     = p0;
exact.B0     = B0;
exact.A      = A;
exact.k      = k;
exact.vA     = vA;
exact.lambda = lambda;
exact.period = period;

exact.uz = @(x,t) A .* sin(k .* (x - vA .* t));
exact.Bz = @(x,t) -sqrt(rho0) .* A .* sin(k .* (x - vA .* t));

params = struct();
params.rho0   = rho0;
params.p0     = p0;
params.B0     = B0;
params.A      = A;
params.mode   = mode;
params.k      = k;
params.vA     = vA;
params.lambda = lambda;
params.period = period;
params.Lx     = Lx;

end
