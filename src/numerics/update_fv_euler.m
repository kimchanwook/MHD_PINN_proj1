function [Unew, stepInfo] = update_fv_euler(U, gridData, gamma, dt, positivity)
% update_fv_euler
%
% Advances the conserved state by one forward-Euler time step:
%
%   Unew = U + dt * rhs
%
% where rhs = dU/dt is computed by compute_rhs.
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
%   stepInfo   - structure with positivity-floor information

if nargin < 5
    positivity = [];
end

rhs  = compute_rhs(U, gridData, gamma);
Unew = U + dt .* rhs;

stepInfo = struct();
stepInfo.floorApplied = false;
stepInfo.floorInfo = [];

if ~isempty(positivity) && isfield(positivity, 'isEnabled') && positivity.isEnabled
    [Unew, floorInfo] = apply_positivity_floors(Unew, gamma, positivity);
    stepInfo.floorApplied = floorInfo.anyFloorApplied;
    stepInfo.floorInfo = floorInfo;
end

end
