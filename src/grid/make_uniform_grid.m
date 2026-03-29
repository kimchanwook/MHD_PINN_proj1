function gridData = make_uniform_grid(xMin, xMax, yMin, yMax, Nx, Ny)
% make_uniform_grid
%
% Creates a uniform 2D Cartesian gridData for a finite-volume solver.

dx = (xMax - xMin) / Nx;
dy = (yMax - yMin) / Ny;

xc = xMin + (0.5 : 1 : Nx - 0.5) * dx;
yc = yMin + (0.5 : 1 : Ny - 0.5)' * dy;

[Xc, Yc] = meshgrid(xc, yc);

xf = xMin + (0 : Nx) * dx;
yf = (yMin + (0 : Ny) * dy)';

gridData = struct();
gridData.xMin = xMin;
gridData.xMax = xMax;
gridData.yMin = yMin;
gridData.yMax = yMax;
gridData.Nx   = Nx;
gridData.Ny   = Ny;
gridData.dx   = dx;
gridData.dy   = dy;
gridData.xc   = xc;
gridData.yc   = yc;
gridData.Xc   = Xc;
gridData.Yc   = Yc;
gridData.xf   = xf;
gridData.yf   = yf;

end
