function residuals = pinn_ideal_mhd_residuals_alfven(net, xCol, yCol, tCol, targets)
% pinn_ideal_mhd_residuals_alfven
%
% Computes reduced residuals for the first executable Alfvén-wave PINN.
%
% IMPORTANT:
%   These residuals are derived from the ideal-MHD equations specialized to
%   the shear-Alfven verification problem used in this project:
%       - fields depend on x and t,
%       - the exact solution is uniform in y,
%       - the dynamical coupling is carried by uz and Bz,
%       - rho, p, and Bx remain at their background values,
%       - ux, uy, and By remain zero.
%
% INPUTS:
%   net     - dlnetwork
%   xCol    - dlarray of x collocation coordinates, size [1, N]
%   yCol    - dlarray of y collocation coordinates, size [1, N]
%   tCol    - dlarray of t collocation coordinates, size [1, N]
%   targets - structure with background values rho0, p0, and B0
%
% OUTPUT:
%   residuals - structure of dlarray residual vectors
%
% REDUCED IDEAL-MHD RESIDUALS:
%   d(uz)/dt - (B0/rho0) * d(Bz)/dx = 0
%   d(Bz)/dt - B0 * d(uz)/dx        = 0
%
% BACKGROUND/Y-UNIFORMITY RESIDUALS:
%   rho - rho0 = 0, p - p0 = 0, Bx - B0 = 0
%   ux = 0, uy = 0, By = 0
%   d()/dy = 0 for all predicted fields
%   div(B) = d(Bx)/dx + d(By)/dy = 0

xyt = [xCol; yCol; tCol];
out = forward(net, xyt);

rho = out(1,:);
ux  = out(2,:);
uy  = out(3,:);
uz  = out(4,:);
p   = out(5,:);
Bx  = out(6,:);
By  = out(7,:);
Bz  = out(8,:);

% Derivatives with respect to x, y, and t.
drho = dlgradient(sum(rho, 'all'), {xCol, yCol, tCol}, 'EnableHigherDerivatives', true);
dux  = dlgradient(sum(ux,  'all'), {xCol, yCol, tCol}, 'EnableHigherDerivatives', true);
duy  = dlgradient(sum(uy,  'all'), {xCol, yCol, tCol}, 'EnableHigherDerivatives', true);
duz  = dlgradient(sum(uz,  'all'), {xCol, yCol, tCol}, 'EnableHigherDerivatives', true);
dp   = dlgradient(sum(p,   'all'), {xCol, yCol, tCol}, 'EnableHigherDerivatives', true);
dBx  = dlgradient(sum(Bx,  'all'), {xCol, yCol, tCol}, 'EnableHigherDerivatives', true);
dBy  = dlgradient(sum(By,  'all'), {xCol, yCol, tCol}, 'EnableHigherDerivatives', true);
dBz  = dlgradient(sum(Bz,  'all'), {xCol, yCol, tCol}, 'EnableHigherDerivatives', true);

rho0 = targets.background.rho0;
p0   = targets.background.p0;
B0   = targets.background.B0;

residuals = struct();

% Reduced shear-Alfven dynamics.
residuals.uzDynamics = duz{3} - (B0 / rho0) .* dBz{1};
residuals.BzDynamics = dBz{3} - B0 .* duz{1};

% Background-state consistency.
residuals.rhoBackground = rho - rho0;
residuals.pBackground   = p   - p0;
residuals.BxBackground  = Bx  - B0;
residuals.uxZero        = ux;
residuals.uyZero        = uy;
residuals.ByZero        = By;

% y-uniformity for the Alfvén benchmark.
residuals.rhoY = drho{2};
residuals.uxY  = dux{2};
residuals.uyY  = duy{2};
residuals.uzY  = duz{2};
residuals.pY   = dp{2};
residuals.BxY  = dBx{2};
residuals.ByY  = dBy{2};
residuals.BzY  = dBz{2};

% Divergence-free magnetic field.
residuals.divB = dBx{1} + dBy{2};

end
