function plot_pulse_fields(grid, V, t, outDir, prefix)
% plot_pulse_fields
%
% Creates field plots for the magnetic-pressure pulse case.
%
% INPUTS:
%   grid   - grid structure
%   V      - primitive state, size [Ny, Nx, 8]
%   t      - current time
%   outDir - output directory
%   prefix - prefix added to output filenames

rho = V(:,:,1);
p   = V(:,:,5);
Bx  = V(:,:,6);
By  = V(:,:,7);
Bz  = V(:,:,8);
ux  = V(:,:,2);
uy  = V(:,:,3);
uz  = V(:,:,4);

Bmag = sqrt(Bx.^2 + By.^2 + Bz.^2);
umag = sqrt(ux.^2 + uy.^2 + uz.^2);

fields = {
    rho, 'Density', 'rho';
    p,   'Gas pressure', 'p';
    Bmag,'Magnetic-field magnitude', 'Bmag';
    umag,'Velocity magnitude', 'umag';
    Bz,  'Out-of-plane magnetic field B_z', 'Bz';
};

for i = 1:size(fields,1)
    data = fields{i,1};
    ttl  = fields{i,2};
    tag  = fields{i,3};

    fig = figure('Visible','off');
    imagesc(grid.xc, grid.yc, data);
    set(gca,'YDir','normal');
    axis equal tight;
    xlabel('x'); ylabel('y');
    title(sprintf('%s at t = %.6f', ttl, t));
    colorbar;
    saveas(fig, fullfile(outDir, sprintf('%s_%s.png', prefix, tag)));
    close(fig);
end

fig = figure('Visible','off');
imagesc(grid.xc, grid.yc, Bmag);
set(gca,'YDir','normal');
hold on;
step = max(1, floor(min(grid.Nx, grid.Ny)/24));
quiver(grid.Xc(1:step:end,1:step:end), grid.Yc(1:step:end,1:step:end), ...
       Bx(1:step:end,1:step:end), By(1:step:end,1:step:end), 'k');
hold off;
axis equal tight;
xlabel('x'); ylabel('y');
title(sprintf('Magnetic magnitude with in-plane field arrows at t = %.6f', t));
colorbar;
saveas(fig, fullfile(outDir, sprintf('%s_Bmag_quiver.png', prefix)));
close(fig);

end
