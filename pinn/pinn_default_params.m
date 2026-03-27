function params = pinn_default_params()
% pinn_default_params
%
% Returns a default parameter structure for the first executable Alfvén-wave
% PINN experiment.
%
% IMPORTANT:
%   This PINN is the first practical training version for the project.
%   It uses the reduced shear-Alfven subsystem derived from ideal MHD for
%   the specific verification case that depends on x and t while remaining
%   uniform in y. It is therefore a case-specific PINN, not yet a fully
%   general 2D ideal-MHD residual engine.

params = struct();

% Network architecture
params.inputSize = 3;          % x, y, t
params.outputSize = 8;         % rho, ux, uy, uz, p, Bx, By, Bz
params.hiddenWidth = 64;
params.numHiddenLayers = 4;

% Sampling counts
params.numCollocation = 3000;
params.numInitial = 1024;
params.numBoundary = 1024;

% Loss weights
params.lambdaIC        = 10.0;
params.lambdaBC        = 2.0;
params.lambdaReducedPDE = 1.0;
params.lambdaYUniform  = 0.5;
params.lambdaBackground = 1.0;
params.lambdaDivB      = 1.0;

% Training options
params.maxEpochs = 1500;
params.learningRate = 1.0e-3;
params.printEvery = 50;
params.gradientDecayFactor = 0.9;
params.squaredGradientDecayFactor = 0.999;
params.miniBatchFraction = 1.0;

% Output and comparison
params.evalTimeFraction = 0.25;
params.makePlots = true;
params.outputDir = fullfile('output', 'pinn_alfven');

end
