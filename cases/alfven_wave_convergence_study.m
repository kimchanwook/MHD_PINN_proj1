function study = alfven_wave_convergence_study()
% alfven_wave_convergence_study
%
% Runs the Alfv'en-wave verification at multiple gridData resolutions and records
% the final L2 error, maximum divergence error, and positivity-floor usage.

clc;

gamma = 5/3;
CFL   = 0.4;
positivity = default_positivity_settings();

NxList = [32, 64, 128, 256];
NyFactor = 4;

xMin = 0.0;
xMax = 1.0;
yMin = 0.0;
yMax = 1.0;

caseParams = struct();
caseParams.rho0 = 1.0;
caseParams.p0   = 1.0;
caseParams.B0   = 1.0;
caseParams.A    = 1.0e-3;
caseParams.mode = 1;

outDir = fullfile('output', 'alfven_wave_convergence');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

nCases = numel(NxList);
errUz = zeros(nCases,1);
errBz = zeros(nCases,1);
maxDivB = zeros(nCases,1);
dxVals = zeros(nCases,1);
measuredSpeed = zeros(nCases,1);
numFloorEvents = zeros(nCases,1);

for iCase = 1:nCases
    Nx = NxList(iCase);
    Ny = max(8, round(Nx / NyFactor));
    gridData = make_uniform_grid(xMin, xMax, yMin, yMax, Nx, Ny);
    [U, exact, params] = init_alfven_wave(gridData, gamma, caseParams);

    V0 = conserved_to_primitive(U, gamma);
    midRow = ceil(gridData.Ny / 2);
    uz_t0 = squeeze(V0(midRow,:,4)).';

    t = 0.0;
    finalTime = params.period;
    step = 0;
    maxSteps = 500000;

    while t < finalTime && step < maxSteps
        dt = compute_time_step(U, gridData, gamma, CFL);
        if t + dt > finalTime
            dt = finalTime - t;
        end
        [U, stepInfo] = update_fv_rk2_plm(U, gridData, gamma, dt, positivity);
        if stepInfo.floorApplied
            numFloorEvents(iCase) = numFloorEvents(iCase) + 1;
        end
        t = t + dt;
        step = step + 1;
    end

    Vfinal = conserved_to_primitive(U, gamma);
    uz_final = squeeze(Vfinal(midRow,:,4)).';
    Bz_final = squeeze(Vfinal(midRow,:,8)).';
    uz_exact = exact.uz(gridData.xc(:), t);
    Bz_exact = exact.Bz(gridData.xc(:), t);
    divB = compute_divB(U, gridData);

    errUz(iCase) = compute_l2_error(uz_final, uz_exact);
    errBz(iCase) = compute_l2_error(Bz_final, Bz_exact);
    maxDivB(iCase) = max(abs(divB(:)));
    dxVals(iCase) = gridData.dx;
    measuredSpeed(iCase) = measure_alfven_phase_speed(gridData.xc, uz_t0, uz_final, t);
end

fig1 = figure('Visible','off');
loglog(dxVals, errUz, 'o-', 'LineWidth', 1.5, 'DisplayName', 'u_z L2 error');
hold on;
loglog(dxVals, errBz, 's-', 'LineWidth', 1.5, 'DisplayName', 'B_z L2 error');
hold off;
set(gca,'XDir','reverse');
xlabel('dx'); ylabel('L2 error');
title('Alfven-wave convergence study');
legend('Location','best'); grid on;
saveas(fig1, fullfile(outDir, 'alfven_convergence.png'));
close(fig1);

fig2 = figure('Visible','off');
loglog(dxVals, maxDivB, 'o-', 'LineWidth', 1.5);
set(gca,'XDir','reverse');
xlabel('dx'); ylabel('max |div B|');
title('Divergence error versus gridData spacing');
grid on;
saveas(fig2, fullfile(outDir, 'alfven_divB_convergence.png'));
close(fig2);

fid = fopen(fullfile(outDir, 'alfven_convergence_summary.txt'), 'w');
fprintf(fid, 'Nx    Ny       dx             L2(uz)         L2(Bz)         max|divB|       measured speed   floor events');
for iCase = 1:nCases
    fprintf(fid, '%-5d %-5d %-14.6e %-14.6e %-14.6e %-14.6e %-16.6f %-8d', ...
        NxList(iCase), max(8, round(NxList(iCase)/NyFactor)), dxVals(iCase), errUz(iCase), errBz(iCase), maxDivB(iCase), measuredSpeed(iCase), numFloorEvents(iCase));
end
fclose(fid);

study = struct();
study.NxList = NxList(:);
study.NyList = arrayfun(@(n) max(8, round(n/NyFactor)), NxList(:));
study.dxVals = dxVals;
study.errUz = errUz;
study.errBz = errBz;
study.maxDivB = maxDivB;
study.measuredSpeed = measuredSpeed;
study.numFloorEvents = numFloorEvents;

end
