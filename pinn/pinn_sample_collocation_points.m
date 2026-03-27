function samples = pinn_sample_collocation_points(grid, tMin, tMax, params)
% pinn_sample_collocation_points
%
% Generates random sample points for the PINN loss terms.
%
% OUTPUT fields:
%   samples.collocation - [Nf, 3] points in interior of space-time domain
%   samples.initial     - [Ni, 3] points on t = tMin plane
%   samples.boundaryL   - [Nb, 3] points on x = xMin plane
%   samples.boundaryR   - [Nb, 3] points on x = xMax plane
%
% NOTES:
%   This first executable PINN uses periodicity only in x for the Alfvén
%   case, because the reference solution depends on x and t while remaining
%   uniform in y.

Nf = params.numCollocation;
Ni = params.numInitial;
Nb = params.numBoundary;

xMin = grid.xMin;
xMax = grid.xMax;
yMin = grid.yMin;
yMax = grid.yMax;

samples = struct();

samples.collocation = [ ...
    xMin + (xMax - xMin) * rand(Nf, 1), ...
    yMin + (yMax - yMin) * rand(Nf, 1), ...
    tMin  + (tMax - tMin) * rand(Nf, 1)  ...
    ];

samples.initial = [ ...
    xMin + (xMax - xMin) * rand(Ni, 1), ...
    yMin + (yMax - yMin) * rand(Ni, 1), ...
    tMin * ones(Ni, 1) ...
    ];

randY = yMin + (yMax - yMin) * rand(Nb, 1);
randT = tMin + (tMax - tMin) * rand(Nb, 1);

samples.boundaryL = [xMin * ones(Nb,1), randY, randT];
samples.boundaryR = [xMax * ones(Nb,1), randY, randT];

end
