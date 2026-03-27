function [U0, params] = init_magnetic_pressure_pulse(grid, gamma, userParams)
% init_magnetic_pressure_pulse
%
% Initializes a localized magnetic-pressure pulse on a uniform plasma
% background. The pulse is imposed in the out-of-plane magnetic-field
% component Bz so that the magnetic energy density is largest near the pulse
% center.
%
% INPUTS:
%   grid       - grid structure from make_uniform_grid
%   gamma      - ratio of specific heats
%   userParams - structure with fields:
%       rho0   - background mass density
%       p0     - background gas pressure
%       Bx0    - background magnetic field in x-direction
%       By0    - background magnetic field in y-direction
%       A      - pulse amplitude added to Bz
%       sigma  - pulse width
%       x0     - pulse center x-coordinate
%       y0     - pulse center y-coordinate
%
% OUTPUTS:
%   U0         - initial conserved state
%   params     - structure containing case parameters

rho0  = userParams.rho0;
p0    = userParams.p0;
Bx0   = userParams.Bx0;
By0   = userParams.By0;
A     = userParams.A;
sigma = userParams.sigma;
x0    = userParams.x0;
y0    = userParams.y0;

Ny = grid.Ny;
Nx = grid.Nx;
nVar = 8;

X = grid.Xc;
Y = grid.Yc;
R2 = (X - x0).^2 + (Y - y0).^2;

V0 = zeros(Ny, Nx, nVar);
V0(:,:,1) = rho0;
V0(:,:,2) = 0.0;
V0(:,:,3) = 0.0;
V0(:,:,4) = 0.0;
V0(:,:,5) = p0;
V0(:,:,6) = Bx0;
V0(:,:,7) = By0;
V0(:,:,8) = A .* exp(-R2 ./ (2 .* sigma.^2));

U0 = primitive_to_conserved(V0, gamma);

params = struct();
params.rho0  = rho0;
params.p0    = p0;
params.Bx0   = Bx0;
params.By0   = By0;
params.A     = A;
params.sigma = sigma;
params.x0    = x0;
params.y0    = y0;

end
