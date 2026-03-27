function grid = make_uniform_grid(xMin, xMax, yMin, yMax, Nx, Ny)
% make_uniform_grid
%
% Creates a uniform 2D Cartesian grid for a finite-volume solver.

dx = (xMax - xMin) / Nx;
dy = (yMax - yMin) / Ny;

xc = xMin + (0.5 : 1 : Nx - 0.5) * dx;
yc = yMin + (0.5 : 1 : Ny - 0.5)' * dy;

[Xc, Yc] = meshgrid(xc, yc);

xf = xMin + (0 : Nx) * dx;
yf = (yMin + (0 : Ny) * dy)';

grid = struct();
grid.xMin = xMin;
grid.xMax = xMax;
grid.yMin = yMin;
grid.yMax = yMax;
grid.Nx   = Nx;
grid.Ny   = Ny;
grid.dx   = dx;
grid.dy   = dy;
grid.xc   = xc;
grid.yc   = yc;
grid.Xc   = Xc;
grid.Yc   = Yc;
grid.xf   = xf;
grid.yf   = yf;

end
