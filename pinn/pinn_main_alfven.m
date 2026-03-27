function pinnResults = pinn_main_alfven()
% pinn_main_alfven
%
% Runs the first executable physics-informed neural network (PINN) training
% loop for the Alfvén-wave verification problem.
%
% IMPORTANT:
%   This is the first practical PINN stage for the project. The residuals are
%   derived from the ideal-MHD equations specialized to the shear-Alfven test
%   case used by the conventional solver. The code is therefore deliberately
%   narrower than a fully general 2D ideal-MHD PINN.

clc;

gamma = 5/3;
params = pinn_default_params();

if ~exist(params.outputDir, 'dir')
    mkdir(params.outputDir);
end

grid = make_uniform_grid(0.0, 1.0, 0.0, 1.0, 128, 32);
caseParams = struct('rho0', 1.0, 'p0', 1.0, 'B0', 1.0, 'A', 1.0e-3, 'mode', 1);
[~, exact, solverParams] = init_alfven_wave(grid, gamma, caseParams);

samples = pinn_sample_collocation_points(grid, 0.0, solverParams.period, params);
targets = pinn_apply_ic_bc(samples, exact);
net = pinn_build_network(params);

initialLosses = pinn_compute_losses(net, samples, targets, params);
fprintf('Initial total loss before training: %.4e\n', initialLosses.total);

trainResults = pinn_train_alfven(net, samples, targets, params);
trainedNet = trainResults.net;

comparison = pinn_compare_with_solver( ...
    trainedNet, grid, exact, params.evalTimeFraction * solverParams.period);

if params.makePlots
    pinn_plot_training_history(trainResults.history, params.outputDir);
    pinn_plot_comparison_alfven(comparison, params.outputDir);
end

summaryFile = fullfile(params.outputDir, 'pinn_alfven_summary.txt');
fid = fopen(summaryFile, 'w');
fprintf(fid, 'Executable PINN summary for the Alfven-wave case\n');
fprintf(fid, '------------------------------------------------\n');
fprintf(fid, 'Network hidden width      = %d\n', params.hiddenWidth);
fprintf(fid, 'Number of hidden layers   = %d\n', params.numHiddenLayers);
fprintf(fid, 'Collocation points        = %d\n', params.numCollocation);
fprintf(fid, 'Initial-condition points  = %d\n', params.numInitial);
fprintf(fid, 'Boundary points           = %d\n', params.numBoundary);
fprintf(fid, 'Max epochs                = %d\n', params.maxEpochs);
fprintf(fid, 'Learning rate             = %.6e\n', params.learningRate);
fprintf(fid, 'Initial total loss        = %.6e\n', initialLosses.total);
fprintf(fid, 'Final total loss          = %.6e\n', trainResults.history.total(end));
fprintf(fid, 'Evaluation time           = %.6f\n', comparison.tEval);
fprintf(fid, 'L2 error in uz            = %.6e\n', comparison.errUzL2);
fprintf(fid, 'L2 error in Bz            = %.6e\n', comparison.errBzL2);
fprintf(fid, 'Max |rho - rho0|          = %.6e\n', comparison.maxAbsRhoDrift);
fprintf(fid, 'Max |p - p0|              = %.6e\n', comparison.maxAbsPDrift);
fprintf(fid, 'Max |Bx - B0|             = %.6e\n', comparison.maxAbsBxDrift);
fprintf(fid, 'Max |uy|                  = %.6e\n', comparison.maxAbsUy);
fprintf(fid, 'Max |By|                  = %.6e\n', comparison.maxAbsBy);
fclose(fid);

fprintf('PINN Alfven training finished.\n');
fprintf('Final total loss: %.4e\n', trainResults.history.total(end));
fprintf('PINN L2 error in u_z: %.4e\n', comparison.errUzL2);
fprintf('PINN L2 error in B_z: %.4e\n', comparison.errBzL2);

pinnResults = struct();
pinnResults.params = params;
pinnResults.grid = grid;
pinnResults.solverParams = solverParams;
pinnResults.samples = samples;
pinnResults.targets = targets;
pinnResults.initialLosses = initialLosses;
pinnResults.trainResults = trainResults;
pinnResults.net = trainedNet;
pinnResults.comparison = comparison;

end
