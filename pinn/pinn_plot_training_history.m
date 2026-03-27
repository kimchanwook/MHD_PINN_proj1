function pinn_plot_training_history(history, outDir)
% pinn_plot_training_history
%
% Saves loss-history plots for the Alfvén-wave PINN training run.

if ~exist(outDir, 'dir')
    mkdir(outDir);
end

fig = figure('Visible', 'off');
semilogy(history.epoch, history.total, 'LineWidth', 1.5, 'DisplayName', 'Total');
hold on;
semilogy(history.epoch, history.initial, 'LineWidth', 1.2, 'DisplayName', 'IC');
semilogy(history.epoch, history.boundary, 'LineWidth', 1.2, 'DisplayName', 'BC');
semilogy(history.epoch, history.reducedPDE, 'LineWidth', 1.2, 'DisplayName', 'Reduced PDE');
semilogy(history.epoch, history.yUniform, 'LineWidth', 1.2, 'DisplayName', 'Y-uniform');
semilogy(history.epoch, history.background, 'LineWidth', 1.2, 'DisplayName', 'Background');
semilogy(history.epoch, history.divB, 'LineWidth', 1.2, 'DisplayName', 'div B');
hold off;
xlabel('Epoch');
ylabel('Loss');
title('PINN training history for the Alfven-wave case');
legend('Location', 'best');
grid on;
saveas(fig, fullfile(outDir, 'pinn_training_history.png'));
close(fig);

end
