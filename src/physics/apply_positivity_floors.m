function [Uout, floorInfo] = apply_positivity_floors(Uin, gamma, positivity)
% apply_positivity_floors
%
% Applies minimum floors to mass density and gas pressure, then rebuilds the
% conserved energy so the primitive and conserved states remain consistent.
%
% INPUTS:
%   Uin        - conserved state array, size [Ny, Nx, 8]
%   gamma      - ratio of specific heats
%   positivity - structure with fields:
%       rhoFloor - minimum allowed mass density
%       pFloor   - minimum allowed gas pressure
%
% OUTPUTS:
%   Uout       - floor-enforced conserved state
%   floorInfo  - structure with counts and extrema before flooring
%
% METHOD:
%   1. Convert Uin to primitive variables.
%   2. Floor rho and p.
%   3. Keep velocities and magnetic field components unchanged.
%   4. Recompute total energy from the floored primitive state.
%
% PHYSICAL INTERPRETATION:
%   This is a numerical safeguard, not a physical model. It is used to avoid
%   negative density or pressure states that can arise from truncation error,
%   overshoots, or insufficient numerical dissipation.

V = conserved_to_primitive(Uin, gamma);

rho = V(:,:,1);
p   = V(:,:,5);

rhoMinBefore = min(rho(:));
pMinBefore   = min(p(:));

rhoFloor = positivity.rhoFloor;
pFloor   = positivity.pFloor;

rhoFloored = max(rho, rhoFloor);
pFloored   = max(p,   pFloor);

rhoFlooredCount = nnz(rho < rhoFloor);
pFlooredCount   = nnz(p   < pFloor);

V(:,:,1) = rhoFloored;
V(:,:,5) = pFloored;

Uout = primitive_to_conserved(V, gamma);

floorInfo = struct();
floorInfo.rhoFloor = rhoFloor;
floorInfo.pFloor = pFloor;
floorInfo.rhoMinBefore = rhoMinBefore;
floorInfo.pMinBefore = pMinBefore;
floorInfo.rhoFlooredCount = rhoFlooredCount;
floorInfo.pFlooredCount = pFlooredCount;
floorInfo.anyFloorApplied = (rhoFlooredCount > 0) || (pFlooredCount > 0);

end
