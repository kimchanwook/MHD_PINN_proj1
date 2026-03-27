function stateInfo = check_physical_state(U, gamma)
% check_physical_state
%
% Computes simple physical-state diagnostics for density and pressure.
%
% INPUTS:
%   U      - conserved state array, size [Ny, Nx, 8]
%   gamma  - ratio of specific heats
%
% OUTPUT:
%   stateInfo - structure with minima and negativity counts

V = conserved_to_primitive(U, gamma);

rho = V(:,:,1);
p   = V(:,:,5);

stateInfo = struct();
stateInfo.rhoMin = min(rho(:));
stateInfo.pMin = min(p(:));
stateInfo.negativeOrZeroRhoCount = nnz(rho <= 0);
stateInfo.negativeOrZeroPCount   = nnz(p   <= 0);

end
