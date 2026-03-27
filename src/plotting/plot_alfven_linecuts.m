function plot_alfven_linecuts(grid, V, exact, t, outDir)
% plot_alfven_linecuts
%
% Creates line-cut plots for the Alfvén wave verification case.

midRow = ceil(grid.Ny / 2);

x = grid.xc(:);
uz_num = squeeze(V(midRow,:,4)).';
Bz_num = squeeze(V(midRow,:,8)).';

uz_exact = exact.uz(x, t);
Bz_exact = exact.Bz(x, t);

fig1 = figure('Visible', 'off');
plot(x, uz_num, 'o-', 'LineWidth', 1.2, 'DisplayName', 'Numerical u_z');
hold on;
plot(x, uz_exact, '--', 'LineWidth', 1.5, 'DisplayName', 'Exact u_z');
hold off;
xlabel('x');
ylabel('u_z');
title(sprintf('Alfven wave line cut: u_z at t = %.6f', t));
legend('Location', 'best');
grid on;
saveas(fig1, fullfile(outDir, 'alfven_linecut_uz.png'));
close(fig1);

fig2 = figure('Visible', 'off');
plot(x, Bz_num, 'o-', 'LineWidth', 1.2, 'DisplayName', 'Numerical B_z');
hold on;
plot(x, Bz_exact, '--', 'LineWidth', 1.5, 'DisplayName', 'Exact B_z');
hold off;
xlabel('x');
ylabel('B_z');
title(sprintf('Alfven wave line cut: B_z at t = %.6f', t));
legend('Location', 'best');
grid on;
saveas(fig2, fullfile(outDir, 'alfven_linecut_Bz.png'));
close(fig2);

end
