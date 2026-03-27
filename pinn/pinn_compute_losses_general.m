function losses = pinn_compute_losses_general(net, samples, targets, params, gamma)
% pinn_compute_losses_general
%
% Computes a broad loss decomposition for a future general ideal-MHD PINN.
% This file is provided as a scaffold for the next project stage.
%
% IMPORTANT:
%   This function does not by itself make the full general PINN practical.
%   It is intentionally a structured bridge between the reduced Alfven PINN
%   and a later broader ideal-MHD PINN.

residuals = pinn_ideal_mhd_residuals_general(net, ...
    samples.collocation.x, samples.collocation.y, samples.collocation.t, gamma);

names = fieldnames(residuals);
pdeLoss = dlarray(0.0);
for k = 1:numel(names)
    r = residuals.(names{k});
    pdeLoss = pdeLoss + mean(r.^2, 'all');
end

icbc = pinn_apply_ic_bc(samples, targets);
baseLosses = pinn_compute_losses(net, samples, icbc, params);

losses = struct();
losses.pdeGeneral = pdeLoss;
losses.icbcBase   = baseLosses.total;
losses.total      = params.lambdaIC * baseLosses.ic + ...
                    params.lambdaBC * baseLosses.bc + ...
                    pdeLoss;
end
