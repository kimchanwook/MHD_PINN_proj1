function [lossValue, gradients, lossStruct] = pinn_model_gradients_alfven(net, samples, targets, params)
% pinn_model_gradients_alfven
%
% Computes the weighted PINN loss and network gradients for the first
% executable Alfvén-wave training loop.
%
% INPUTS:
%   net     - dlnetwork
%   samples - sample structure from pinn_sample_collocation_points
%   targets - target structure from pinn_apply_ic_bc
%   params  - PINN parameter structure
%
% OUTPUTS:
%   lossValue - scalar dlarray total loss
%   gradients - dlgradient of total loss with respect to learnables
%   lossStruct - structure of individual scalar losses

% -------------------------------------------------------------------------
% Initial-condition loss
% -------------------------------------------------------------------------
Xi = dlarray(samples.initial(:,1).', 'CB');
Yi = dlarray(samples.initial(:,2).', 'CB');
Ti = dlarray(samples.initial(:,3).', 'CB');

outIC = forward(net, [Xi; Yi; Ti]);

lossIC = mean((outIC(1,:) - targets.initial.rho.').^2, 'all') + ...
         mean((outIC(2,:) - targets.initial.ux.').^2,  'all') + ...
         mean((outIC(3,:) - targets.initial.uy.').^2,  'all') + ...
         mean((outIC(4,:) - targets.initial.uz.').^2,  'all') + ...
         mean((outIC(5,:) - targets.initial.p.').^2,   'all') + ...
         mean((outIC(6,:) - targets.initial.Bx.').^2,  'all') + ...
         mean((outIC(7,:) - targets.initial.By.').^2,  'all') + ...
         mean((outIC(8,:) - targets.initial.Bz.').^2,  'all');

% -------------------------------------------------------------------------
% Periodic boundary loss in x
% -------------------------------------------------------------------------
XL = dlarray(samples.boundaryL(:,1).', 'CB');
YL = dlarray(samples.boundaryL(:,2).', 'CB');
TL = dlarray(samples.boundaryL(:,3).', 'CB');

XR = dlarray(samples.boundaryR(:,1).', 'CB');
YR = dlarray(samples.boundaryR(:,2).', 'CB');
TR = dlarray(samples.boundaryR(:,3).', 'CB');

outL = forward(net, [XL; YL; TL]);
outR = forward(net, [XR; YR; TR]);
lossBC = mean((outL - outR).^2, 'all');

% -------------------------------------------------------------------------
% Reduced ideal-MHD residual loss on collocation points
% -------------------------------------------------------------------------
Xc = dlarray(samples.collocation(:,1).', 'CB');
Yc = dlarray(samples.collocation(:,2).', 'CB');
Tc = dlarray(samples.collocation(:,3).', 'CB');

residuals = pinn_ideal_mhd_residuals_alfven(net, Xc, Yc, Tc, targets);

lossReducedPDE = mean(residuals.uzDynamics.^2, 'all') + ...
                 mean(residuals.BzDynamics.^2, 'all');

lossYUniform = mean(residuals.rhoY.^2, 'all') + ...
               mean(residuals.uxY.^2,  'all') + ...
               mean(residuals.uyY.^2,  'all') + ...
               mean(residuals.uzY.^2,  'all') + ...
               mean(residuals.pY.^2,   'all') + ...
               mean(residuals.BxY.^2,  'all') + ...
               mean(residuals.ByY.^2,  'all') + ...
               mean(residuals.BzY.^2,  'all');

lossBackground = mean(residuals.rhoBackground.^2, 'all') + ...
                 mean(residuals.pBackground.^2,   'all') + ...
                 mean(residuals.BxBackground.^2,  'all') + ...
                 mean(residuals.uxZero.^2,        'all') + ...
                 mean(residuals.uyZero.^2,        'all') + ...
                 mean(residuals.ByZero.^2,        'all');

lossDivB = mean(residuals.divB.^2, 'all');

lossValue = params.lambdaIC         * lossIC + ...
            params.lambdaBC         * lossBC + ...
            params.lambdaReducedPDE * lossReducedPDE + ...
            params.lambdaYUniform   * lossYUniform + ...
            params.lambdaBackground * lossBackground + ...
            params.lambdaDivB       * lossDivB;

lossStruct = struct();
lossStruct.initial = gather(extractdata(lossIC));
lossStruct.boundary = gather(extractdata(lossBC));
lossStruct.reducedPDE = gather(extractdata(lossReducedPDE));
lossStruct.yUniform = gather(extractdata(lossYUniform));
lossStruct.background = gather(extractdata(lossBackground));
lossStruct.divB = gather(extractdata(lossDivB));
lossStruct.total = gather(extractdata(lossValue));

gradients = dlgradient(lossValue, net.Learnables);

end
