function pinn_plot_comparison_alfven(comparison, outDir)
% pinn_plot_comparison_alfven
%
% Saves comparison plots between the PINN prediction and the analytical
% Alfvén-wave reference.

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

fig1 = figure('Visible', 'off');
plot(comparison.x, comparison.uzPred, 'o-', 'LineWidth', 1.2, 'DisplayName', 'PINN u_z');
hold on;
plot(comparison.x, comparison.uzExact, '--', 'LineWidth', 1.5, 'DisplayName', 'Exact u_z');
hold off;
xlabel('x');
ylabel('u_z');
title(sprintf('PINN vs exact: u_z at t = %.6f', comparison.tEval));
legend('Location', 'best');
grid on;
saveas(fig1, fullfile(outDir, 'pinn_compare_uz.png'));
close(fig1);

fig2 = figure('Visible', 'off');
plot(comparison.x, comparison.BzPred, 'o-', 'LineWidth', 1.2, 'DisplayName', 'PINN B_z');
hold on;
plot(comparison.x, comparison.BzExact, '--', 'LineWidth', 1.5, 'DisplayName', 'Exact B_z');
hold off;
xlabel('x');
ylabel('B_z');
title(sprintf('PINN vs exact: B_z at t = %.6f', comparison.tEval));
legend('Location', 'best');
grid on;
saveas(fig2, fullfile(outDir, 'pinn_compare_Bz.png'));
close(fig2);

end
