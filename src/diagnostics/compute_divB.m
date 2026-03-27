function [divB, maxAbsDivB, rmsDivB] = compute_divB(U, grid)
% compute_divB
%
% Computes the magnetic-field divergence diagnostic on the physical grid.
%
% INPUTS:
%   U      - conserved state on physical domain, size [Ny, Nx, 8]
%   grid   - grid structure from make_uniform_grid
%
% OUTPUTS:
%   divB       - discrete divergence of magnetic field on cell centers
%   maxAbsDivB - maximum absolute divergence over the domain
%   rmsDivB    - root-mean-square divergence over the domain
%
% DEFINITION:
%   divB = dBx/dx + dBy/dy
%
% DISCRETIZATION:
%   Uses second-order central differences with periodic wrapping:
%
%     dBx/dx ~ [Bx(i+1) - Bx(i-1)] / (2*dx)
%     dBy/dy ~ [By(j+1) - By(j-1)] / (2*dy)
%
% PHYSICAL COMMENT:
%   Ideal MHD requires nabla · B = 0. Numerically, this is rarely exactly
%   satisfied unless a divergence-control strategy is used. This routine is
%   a diagnostic that tells us how large the magnetic-monopole error has
%   become.

Bx = U(:,:,6);
By = U(:,:,7);

dBx_dx = (circshift(Bx, [0, -1]) - circshift(Bx, [0, 1])) ./ (2 .* grid.dx);
dBy_dy = (circshift(By, [-1, 0]) - circshift(By, [1, 0])) ./ (2 .* grid.dy);

divB = dBx_dx + dBy_dy;
maxAbsDivB = max(abs(divB(:)));
rmsDivB = sqrt(mean(divB(:).^2));

end
