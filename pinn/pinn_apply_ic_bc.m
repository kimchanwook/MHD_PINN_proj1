function targets = pinn_apply_ic_bc(samples, exact)
% pinn_apply_ic_bc
%
% Builds target values for the initial-condition and periodic-boundary parts
% of the Alfvén-wave PINN problem.
%
% INPUTS:
%   samples - structure from pinn_sample_collocation_points
%   exact   - exact solution structure from init_alfven_wave
%
% OUTPUT:
%   targets - structure with target values for loss construction

Xi = samples.initial(:,1);
Ti = samples.initial(:,3);

targets = struct();

% Initial-condition targets
Ninit = size(samples.initial, 1);
targets.initial = struct();
targets.initial.rho = exact.rho0 * ones(Ninit, 1);
targets.initial.ux  = zeros(Ninit, 1);
targets.initial.uy  = zeros(Ninit, 1);
targets.initial.uz  = exact.uz(Xi, Ti);
targets.initial.p   = exact.p0 * ones(Ninit, 1);
targets.initial.Bx  = exact.B0 * ones(Ninit, 1);
targets.initial.By  = zeros(Ninit, 1);
targets.initial.Bz  = exact.Bz(Xi, Ti);

% Periodic boundary data are matched left-to-right, so only point counts are
% needed here. The actual values come from the network predictions.
Nbc = size(samples.boundaryL, 1);
targets.boundary = struct();
targets.boundary.numPoints = Nbc;

% Background constants used by the reduced Alfvén residuals.
targets.background = struct();
targets.background.rho0 = exact.rho0;
targets.background.p0   = exact.p0;
targets.background.B0   = exact.B0;

end
