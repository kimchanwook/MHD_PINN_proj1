function dt = compute_time_step(U, gridData, gamma, CFL)
% compute_time_step
%
% Computes the explicit time step for ideal MHD using a CFL condition.

rho = U(:,:,1);
mx  = U(:,:,2);
my  = U(:,:,3);

ux = mx ./ rho;
uy = my ./ rho;

cfx = fast_magnetosonic_speed_x(U, gamma);
cfy = fast_magnetosonic_speed_y(U, gamma);

ax = abs(ux) + cfx;
ay = abs(uy) + cfy;

max_ax = max(ax(:));
max_ay = max(ay(:));

dt_x = gridData.dx / max_ax;
dt_y = gridData.dy / max_ay;

dt = CFL * min(dt_x, dt_y);

end
