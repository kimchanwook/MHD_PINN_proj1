function residuals = pinn_ideal_mhd_residuals_general(net, xCol, yCol, tCol, gamma)
% pinn_ideal_mhd_residuals_general
%
% Computes a broad ideal-MHD residual set for a 2D (2.5D) physics-informed
% neural network. This file is meant to be the next step beyond the reduced
% shear-Alfven residuals used for the first executable PINN benchmark.
%
% IMPORTANT:
%   This function is intentionally written as a residual engine, not yet as a
%   complete training workflow. A fully general ideal-MHD PINN is much more
%   demanding than the reduced Alfven case because the equations are strongly
%   nonlinear, coupled, and sensitive to scaling and loss weighting.
%
% INPUTS:
%   net   - dlnetwork
%   xCol  - dlarray of x coordinates, size [1, N]
%   yCol  - dlarray of y coordinates, size [1, N]
%   tCol  - dlarray of t coordinates, size [1, N]
%   gamma - ratio of specific heats
%
% OUTPUT:
%   residuals - structure containing residual vectors for the primitive-form
%               ideal-MHD equations in 2D (2.5D)
%
% MODEL OUTPUT ORDER:
%   rho, ux, uy, uz, p, Bx, By, Bz
%
% EQUATIONS USED:
%   1. Continuity:
%      d(rho)/dt + d(rho*ux)/dx + d(rho*uy)/dy = 0
%
%   2. Momentum in x:
%      d(ux)/dt + ux*d(ux)/dx + uy*d(ux)/dy
%        + (1/rho)*d(ptot)/dx
%        - (1/rho)*( Bx*d(Bx)/dx + By*d(Bx)/dy + Bz*d(Bx? no) )
%
%      In practice it is cleaner to evaluate the conservative momentum-flux
%      form and divide by none; however, for PINN implementation with network
%      outputs in primitive variables, a primitive-form residual is often more
%      straightforward. We therefore use the material-derivative form with
%      magnetic-stress terms derived componentwise below.
%
%   3. Momentum in y
%   4. Momentum in z
%   5. Pressure (adiabatic form):
%      d(p)/dt + ux*d(p)/dx + uy*d(p)/dy + gamma*p*(d(ux)/dx + d(uy)/dy) = 0
%   6. Induction for Bx:
%      d(Bx)/dt + d(uy*Bx - ux*By)/dy = 0
%   7. Induction for By:
%      d(By)/dt + d(ux*By - uy*Bx)/dx = 0
%   8. Induction for Bz:
%      d(Bz)/dt + d(ux*Bz - uz*Bx)/dx + d(uy*Bz - uz*By)/dy = 0
%   9. Divergence-free condition:
%      d(Bx)/dx + d(By)/dy = 0
%
% NOTES:
%   - The pressure equation is the ideal adiabatic closure in primitive form.
%   - The momentum equations below are written using the Lorentz-force form
%     J x B = (curl B) x B, expanded for 2D variation in x and y.
%   - Because Bz can vary with x and y, it contributes to transverse tension.

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

% First derivatives.
drho = dlgradient(sum(rho, 'all'), {xCol, yCol, tCol}, 'EnableHigherDerivatives', true);
dux  = dlgradient(sum(ux,  'all'), {xCol, yCol, tCol}, 'EnableHigherDerivatives', true);
duy  = dlgradient(sum(uy,  'all'), {xCol, yCol, tCol}, 'EnableHigherDerivatives', true);
duz  = dlgradient(sum(uz,  'all'), {xCol, yCol, tCol}, 'EnableHigherDerivatives', true);
dp   = dlgradient(sum(p,   'all'), {xCol, yCol, tCol}, 'EnableHigherDerivatives', true);
dBx  = dlgradient(sum(Bx,  'all'), {xCol, yCol, tCol}, 'EnableHigherDerivatives', true);
dBy  = dlgradient(sum(By,  'all'), {xCol, yCol, tCol}, 'EnableHigherDerivatives', true);
dBz  = dlgradient(sum(Bz,  'all'), {xCol, yCol, tCol}, 'EnableHigherDerivatives', true);

rho_x = drho{1}; rho_y = drho{2}; rho_t = drho{3};
ux_x  = dux{1};  ux_y  = dux{2};  ux_t  = dux{3};
uy_x  = duy{1};  uy_y  = duy{2};  uy_t  = duy{3};
uz_x  = duz{1};  uz_y  = duz{2};  uz_t  = duz{3};
p_x   = dp{1};   p_y   = dp{2};   p_t   = dp{3};
Bx_x  = dBx{1};  Bx_y  = dBx{2};  Bx_t  = dBx{3};
By_x  = dBy{1};  By_y  = dBy{2};  By_t  = dBy{3};
Bz_x  = dBz{1};  Bz_y  = dBz{2};  Bz_t  = dBz{3};

% Derived quantities.
div_u = ux_x + uy_y;
div_B = Bx_x + By_y;

% In 2D (with variation only in x and y), the current-density components in
% normalized units are:
%   Jx =  d(Bz)/dy
%   Jy = -d(Bz)/dx
%   Jz =  d(By)/dx - d(Bx)/dy
Jx = Bz_y;
Jy = -Bz_x;
Jz = By_x - Bx_y;

% Lorentz force J x B.
FxL = Jy .* Bz - Jz .* By;
FyL = Jz .* Bx - Jx .* Bz;
FzL = Jx .* By - Jy .* Bx;

residuals = struct();

% Continuity equation.
residuals.continuity = rho_t + (rho_x .* ux + rho .* ux_x) + (rho_y .* uy + rho .* uy_y);

% Momentum equations in primitive form.
residuals.momentumX = rho .* (ux_t + ux .* ux_x + uy .* ux_y) + p_x - FxL;
residuals.momentumY = rho .* (uy_t + ux .* uy_x + uy .* uy_y) + p_y - FyL;
residuals.momentumZ = rho .* (uz_t + ux .* uz_x + uy .* uz_y)       - FzL;

% Adiabatic pressure equation.
residuals.pressure = p_t + ux .* p_x + uy .* p_y + gamma .* p .* div_u;

% Induction equations in conservative flux form.
fluxBx_y = uy .* Bx - ux .* By;
fluxBy_x = ux .* By - uy .* Bx;
fluxBz_x = ux .* Bz - uz .* Bx;
fluxBz_y = uy .* Bz - uz .* By;

dFluxBx_y = dlgradient(sum(fluxBx_y, 'all'), yCol, 'EnableHigherDerivatives', true);
dFluxBy_x = dlgradient(sum(fluxBy_x, 'all'), xCol, 'EnableHigherDerivatives', true);
dFluxBz_x = dlgradient(sum(fluxBz_x, 'all'), xCol, 'EnableHigherDerivatives', true);
dFluxBz_y = dlgradient(sum(fluxBz_y, 'all'), yCol, 'EnableHigherDerivatives', true);

residuals.inductionBx = Bx_t + dFluxBx_y;
residuals.inductionBy = By_t + dFluxBy_x;
residuals.inductionBz = Bz_t + dFluxBz_x + dFluxBz_y;

% Divergence-free constraint.
residuals.divB = div_B;

% Auxiliary residuals that are often useful during development.
residuals.densityPositive  = -rho;
residuals.pressurePositive = -p;

end
