function Gnum = rusanov_flux_y(UL, UR, gamma)
% rusanov_flux_y
%
% Computes the Rusanov numerical flux in the y-direction for ideal MHD.

GyL = flux_y(UL, gamma);
GyR = flux_y(UR, gamma);

rhoL = UL(:,:,1);
myL  = UL(:,:,3);
uyL  = myL ./ rhoL;
cfyL = fast_magnetosonic_speed_y(UL, gamma);
ayL  = abs(uyL) + cfyL;

rhoR = UR(:,:,1);
myR  = UR(:,:,3);
uyR  = myR ./ rhoR;
cfyR = fast_magnetosonic_speed_y(UR, gamma);
ayR  = abs(uyR) + cfyR;

amax = max(ayL, ayR);

Gnum = 0.5 .* (GyL + GyR) - 0.5 .* amax .* (UR - UL);

end
