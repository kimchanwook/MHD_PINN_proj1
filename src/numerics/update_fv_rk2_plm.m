function [Unew, stepInfo] = update_fv_rk2_plm(U, gridData, gamma, dt, positivity)
% update_fv_rk2_plm
%
% Advances the conserved state by one second-order Runge-Kutta step using
% piecewise-linear reconstruction with a minmod limiter.
%
% INPUTS:
%   U          - conserved state on physical domain, size [Ny, Nx, 8]
%   gridData       - gridData structure
%   gamma      - ratio of specific heats
%   dt         - time-step size
%   positivity - optional positivity-floor structure; omit or pass [] to
%                disable floor enforcement for this update
%
% OUTPUTS:
%   Unew       - updated conserved state
%   stepInfo   - structure describing any floor applications

if nargin < 5
    positivity = [];
end

rhs1 = compute_rhs_plm(U, gridData, gamma);
U1 = U + dt .* rhs1;

stage1FloorInfo = [];
if ~isempty(positivity) && isfield(positivity, 'isEnabled') && positivity.isEnabled
    [U1, stage1FloorInfo] = apply_positivity_floors(U1, gamma, positivity);
end

rhs2 = compute_rhs_plm(U1, gridData, gamma);
Unew = 0.5 .* (U + U1 + dt .* rhs2);

finalFloorInfo = [];
if ~isempty(positivity) && isfield(positivity, 'isEnabled') && positivity.isEnabled
    [Unew, finalFloorInfo] = apply_positivity_floors(Unew, gamma, positivity);
end

stepInfo = struct();
stepInfo.stage1FloorInfo = stage1FloorInfo;
stepInfo.finalFloorInfo = finalFloorInfo;
stepInfo.floorApplied = false;

if ~isempty(stage1FloorInfo)
    stepInfo.floorApplied = stepInfo.floorApplied || stage1FloorInfo.anyFloorApplied;
end
if ~isempty(finalFloorInfo)
    stepInfo.floorApplied = stepInfo.floorApplied || finalFloorInfo.anyFloorApplied;
end

end
