function Fnum = rusanov_flux_x(UL, UR, gamma)
% rusanov_flux_x
%
% Computes the Rusanov numerical flux in the x-direction for ideal MHD.

FxL = flux_x(UL, gamma);
FxR = flux_x(UR, gamma);

rhoL = UL(:,:,1);
mxL  = UL(:,:,2);
uxL  = mxL ./ rhoL;
cfxL = fast_magnetosonic_speed_x(UL, gamma);
axL  = abs(uxL) + cfxL;

rhoR = UR(:,:,1);
mxR  = UR(:,:,2);
uxR  = mxR ./ rhoR;
cfxR = fast_magnetosonic_speed_x(UR, gamma);
axR  = abs(uxR) + cfxR;

amax = max(axL, axR);

Fnum = 0.5 .* (FxL + FxR) - 0.5 .* amax .* (UR - UL);

end
