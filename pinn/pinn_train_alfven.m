function trainResults = pinn_train_alfven(net, samples, targets, params)
% pinn_train_alfven
%
% Trains the Alfvén-wave PINN with Adam using the reduced shear-Alfven
% residual system derived from ideal MHD.
%
% INPUTS:
%   net     - dlnetwork
%   samples - sample structure
%   targets - target structure
%   params  - PINN parameter structure
%
% OUTPUT:
%   trainResults - structure containing the trained network and loss history

trailingAvg = [];
trailingAvgSq = [];

history.epoch = zeros(params.maxEpochs, 1);
history.total = zeros(params.maxEpochs, 1);
history.initial = zeros(params.maxEpochs, 1);
history.boundary = zeros(params.maxEpochs, 1);
history.reducedPDE = zeros(params.maxEpochs, 1);
history.yUniform = zeros(params.maxEpochs, 1);
history.background = zeros(params.maxEpochs, 1);
history.divB = zeros(params.maxEpochs, 1);

for epoch = 1:params.maxEpochs
    [lossValue, gradients, lossStruct] = dlfeval(@pinn_model_gradients_alfven, net, samples, targets, params);

    [net, trailingAvg, trailingAvgSq] = adamupdate( ...
        net, gradients, trailingAvg, trailingAvgSq, epoch, ...
        params.learningRate, params.gradientDecayFactor, params.squaredGradientDecayFactor);

    history.epoch(epoch) = epoch;
    history.total(epoch) = lossStruct.total;
    history.initial(epoch) = lossStruct.initial;
    history.boundary(epoch) = lossStruct.boundary;
    history.reducedPDE(epoch) = lossStruct.reducedPDE;
    history.yUniform(epoch) = lossStruct.yUniform;
    history.background(epoch) = lossStruct.background;
    history.divB(epoch) = lossStruct.divB;

    if mod(epoch, params.printEvery) == 0 || epoch == 1 || epoch == params.maxEpochs
        fprintf(['Epoch %5d | total %.4e | IC %.4e | BC %.4e | PDE %.4e | ' ...
                 'Y %.4e | BG %.4e | divB %.4e\n'], ...
                 epoch, lossStruct.total, lossStruct.initial, lossStruct.boundary, ...
                 lossStruct.reducedPDE, lossStruct.yUniform, ...
                 lossStruct.background, lossStruct.divB);
    end
end

trainResults = struct();
trainResults.net = net;
trainResults.history = history;

end
