function [ULx, URx, ULy, URy] = reconstruct_plm_minmod(U)
% reconstruct_plm_minmod
%
% Piecewise-linear reconstruction with a minmod slope limiter and one layer
% of periodic ghost cells.
%
% INPUT:
%   U      - conserved state on physical domain, size [Ny, Nx, 8]
%
% OUTPUTS:
%   ULx    - left states at x-interfaces,  size [Ny, Nx+1, 8]
%   URx    - right states at x-interfaces, size [Ny, Nx+1, 8]
%   ULy    - lower states at y-interfaces, size [Ny+1, Nx, 8]
%   URy    - upper states at y-interfaces, size [Ny+1, Nx, 8]
%
% METHOD:
%   1. Add one layer of periodic ghost cells.
%   2. Compute limited slopes in x and y using neighboring cell averages.
%   3. Extrapolate each cell-centered value to its faces.
%
% NOTE:
%   This reconstruction is performed directly on conserved variables. That is
%   common in simple finite-volume prototypes, though more advanced solvers
%   often reconstruct primitive or characteristic variables.

[Ny, Nx, nVar] = size(U);
Ubc = apply_periodic_bc(U);

% Limited slopes for each physical cell.
dUx_minus = Ubc(2:Ny+1, 2:Nx+1, :) - Ubc(2:Ny+1, 1:Nx,   :);
dUx_plus  = Ubc(2:Ny+1, 3:Nx+2, :) - Ubc(2:Ny+1, 2:Nx+1, :);
slopeX = minmod(dUx_minus, dUx_plus);

dUy_minus = Ubc(2:Ny+1, 2:Nx+1, :) - Ubc(1:Ny,   2:Nx+1, :);
dUy_plus  = Ubc(3:Ny+2, 2:Nx+1, :) - Ubc(2:Ny+1, 2:Nx+1, :);
slopeY = minmod(dUy_minus, dUy_plus);

% Reconstruct x-interface states.
ULx = zeros(Ny, Nx+1, nVar);
URx = zeros(Ny, Nx+1, nVar);

% Interface 1: between last physical cell and first physical cell (periodic)
ULx(:,1,:)    = Ubc(2:Ny+1, Nx+1, :) + 0.5 .* slopeX(:,Nx,:);
URx(:,1,:)    = Ubc(2:Ny+1, 2,    :) - 0.5 .* slopeX(:,1,:);

% Interior interfaces i+1/2 for i = 1,...,Nx-1
ULx(:,2:Nx,:) = Ubc(2:Ny+1, 2:Nx,   :) + 0.5 .* slopeX(:,1:Nx-1,:);
URx(:,2:Nx,:) = Ubc(2:Ny+1, 3:Nx+1, :) - 0.5 .* slopeX(:,2:Nx,:);

% Interface Nx+1: periodic copy of interface 1 for convenience
ULx(:,Nx+1,:) = ULx(:,1,:);
URx(:,Nx+1,:) = URx(:,1,:);

% Reconstruct y-interface states.
ULy = zeros(Ny+1, Nx, nVar);
URy = zeros(Ny+1, Nx, nVar);

% Interface 1: between last physical row and first physical row (periodic)
ULy(1,:,:)    = Ubc(Ny+1, 2:Nx+1, :) + 0.5 .* slopeY(Ny,:,:);
URy(1,:,:)    = Ubc(2,    2:Nx+1, :) - 0.5 .* slopeY(1,:,:);

% Interior interfaces j+1/2 for j = 1,...,Ny-1
ULy(2:Ny,:,:) = Ubc(2:Ny,   2:Nx+1, :) + 0.5 .* slopeY(1:Ny-1,:,:);
URy(2:Ny,:,:) = Ubc(3:Ny+1, 2:Nx+1, :) - 0.5 .* slopeY(2:Ny,:,:);

% Interface Ny+1: periodic copy of interface 1
ULy(Ny+1,:,:) = ULy(1,:,:);
URy(Ny+1,:,:) = URy(1,:,:);

end
