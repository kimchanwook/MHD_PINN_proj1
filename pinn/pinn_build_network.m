function net = pinn_build_network(params)
% pinn_build_network
%
% Builds a fully connected neural network for the Alfvén-wave PINN.
%
% INPUT:
%   params - structure from pinn_default_params
%
% OUTPUT:
%   net    - dlnetwork object

layers = [ ...
    featureInputLayer(params.inputSize, 'Name', 'input', 'Normalization', 'none')
    fullyConnectedLayer(params.hiddenWidth, 'Name', 'fc1')
    tanhLayer('Name', 'tanh1')
    ];

for k = 2:params.numHiddenLayers
    layers = [layers; ...
        fullyConnectedLayer(params.hiddenWidth, 'Name', sprintf('fc%d', k))
        tanhLayer('Name', sprintf('tanh%d', k))
        ];
end

layers = [layers; ...
    fullyConnectedLayer(params.outputSize, 'Name', 'output')
    ];

lgraph = layerGraph(layers);
net = dlnetwork(lgraph);

end
