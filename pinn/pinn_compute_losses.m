function losses = pinn_compute_losses(net, samples, targets, params)
% pinn_compute_losses
%
% Evaluates the scalar loss components without updating the network.
%
% INPUTS:
%   net     - dlnetwork
%   samples - sample structure
%   targets - target structure
%   params  - PINN parameter structure
%
% OUTPUT:
%   losses  - structure with the same scalar fields reported during training

[~, ~, losses] = dlfeval(@pinn_model_gradients_alfven, net, samples, targets, params);

end
