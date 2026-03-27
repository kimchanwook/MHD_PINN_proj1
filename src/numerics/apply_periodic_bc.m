function Ubc = apply_periodic_bc(U)
% apply_periodic_bc
%
% Adds one layer of periodic ghost cells around a 2D conserved state array.

[Ny, Nx, nVar] = size(U);

Ubc = zeros(Ny + 2, Nx + 2, nVar);

Ubc(2:Ny+1, 2:Nx+1, :) = U;

Ubc(2:Ny+1, 1,      :) = U(:, Nx, :);
Ubc(2:Ny+1, Nx+2,   :) = U(:, 1,  :);

Ubc(1,      2:Nx+1, :) = U(Ny, :, :);
Ubc(Ny+2,   2:Nx+1, :) = U(1,  :, :);

Ubc(1,    1,    :) = U(Ny, Nx, :);
Ubc(1,    Nx+2, :) = U(Ny, 1,  :);
Ubc(Ny+2, 1,    :) = U(1,  Nx, :);
Ubc(Ny+2, Nx+2, :) = U(1,  1,  :);

end
