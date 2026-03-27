function results = alfven_wave_case()
% alfven_wave_case
%
% Runs the Alfv'en-wave verification case using second-order Runge-Kutta time
% stepping and piecewise-linear reconstruction with a minmod limiter.
%
% This stage also enables optional positivity floors for density and gas
% pressure so the solver has a basic robustness safeguard.
%
% OUTPUT:
%   results - structure containing diagnostics and final solution

clc;

gamma = 5/3;
CFL   = 0.4;

Nx = 128;
Ny = 32;

xMin = 0.0;
xMax = 1.0;
yMin = 0.0;
yMax = 1.0;

grid = make_uniform_grid(xMin, xMax, yMin, yMax, Nx, Ny);
positivity = default_positivity_settings();

caseParams = struct();
caseParams.rho0 = 1.0;
caseParams.p0   = 1.0;
caseParams.B0   = 1.0;
caseParams.A    = 1.0e-3;
caseParams.mode = 1;

outDir = fullfile('output', 'alfven_wave');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

[U, exact, params] = init_alfven_wave(grid, gamma, caseParams);
V0 = conserved_to_primitive(U, gamma);
midRow = ceil(grid.Ny / 2);
uz_t0 = squeeze(V0(midRow,:,4)).';

finalTime = params.period;
t = 0.0;
step = 0;
maxSteps = 200000;

numFloorEvents = 0;
numRhoFloorCells = 0;
numPFloorCells = 0;

timeHistory = [];
dtHistory = [];
errUzHistory = [];
errBzHistory = [];
divBmaxHistory = [];

while t < finalTime && step < maxSteps
    dt = compute_time_step(U, grid, gamma, CFL);
    if t + dt > finalTime
        dt = finalTime - t;
    end

    [U, stepInfo] = update_fv_rk2_plm(U, grid, gamma, dt, positivity);

    if stepInfo.floorApplied
        numFloorEvents = numFloorEvents + 1;
    end
    if ~isempty(stepInfo.stage1FloorInfo)
        numRhoFloorCells = numRhoFloorCells + stepInfo.stage1FloorInfo.rhoFlooredCount;
        numPFloorCells   = numPFloorCells   + stepInfo.stage1FloorInfo.pFlooredCount;
    end
    if ~isempty(stepInfo.finalFloorInfo)
        numRhoFloorCells = numRhoFloorCells + stepInfo.finalFloorInfo.rhoFlooredCount;
        numPFloorCells   = numPFloorCells   + stepInfo.finalFloorInfo.pFlooredCount;
    end

    t = t + dt;
    step = step + 1;

    V = conserved_to_primitive(U, gamma);
    uz_num = squeeze(V(midRow,:,4)).';
    Bz_num = squeeze(V(midRow,:,8)).';

    uz_exact = exact.uz(grid.xc(:), t);
    Bz_exact = exact.Bz(grid.xc(:), t);

    errUz = compute_l2_error(uz_num, uz_exact);
    errBz = compute_l2_error(Bz_num, Bz_exact);
    divB = compute_divB(U, grid);

    timeHistory(end+1,1) = t; %#ok<AGROW>
    dtHistory(end+1,1) = dt; %#ok<AGROW>
    errUzHistory(end+1,1) = errUz; %#ok<AGROW>
    errBzHistory(end+1,1) = errBz; %#ok<AGROW>
    divBmaxHistory(end+1,1) = max(abs(divB(:))); %#ok<AGROW>
end

if step >= maxSteps
    warning('Maximum step count reached before final time.');
end

Vfinal = conserved_to_primitive(U, gamma);
uz_final = squeeze(Vfinal(midRow,:,4)).';
Bz_final = squeeze(Vfinal(midRow,:,8)).';
uz_exact_final = exact.uz(grid.xc(:), t);
Bz_exact_final = exact.Bz(grid.xc(:), t);

finalErrUz = compute_l2_error(uz_final, uz_exact_final);
finalErrBz = compute_l2_error(Bz_final, Bz_exact_final);
measuredSpeed = measure_alfven_phase_speed(grid.xc, uz_t0, uz_final, t);
finalDivB = compute_divB(U, grid);
maxAbsDivB = max(abs(finalDivB(:)));
stateInfo = check_physical_state(U, gamma);

plot_alfven_linecuts(grid, Vfinal, exact, t, outDir);

figErr = figure('Visible', 'off');
plot(timeHistory, errUzHistory, 'LineWidth', 1.5, 'DisplayName', 'L2 error in u_z');
hold on;
plot(timeHistory, errBzHistory, 'LineWidth', 1.5, 'DisplayName', 'L2 error in B_z');
hold off;
xlabel('time'); ylabel('L2 error');
title('Alfven-wave error history');
legend('Location', 'best'); grid on;
saveas(figErr, fullfile(outDir, 'alfven_error_history.png'));
close(figErr);

figDiv = figure('Visible', 'off');
plot(timeHistory, divBmaxHistory, 'LineWidth', 1.5);
xlabel('time'); ylabel('max |div B|');
title('Alfven-wave divergence diagnostic');
grid on;
saveas(figDiv, fullfile(outDir, 'alfven_divB_history.png'));
close(figDiv);

summaryFile = fullfile(outDir, 'alfven_summary.txt');
fid = fopen(summaryFile, 'w');
fprintf(fid, 'Alfven wave verification summary
');
fprintf(fid, '--------------------------------
');
fprintf(fid, 'Spatial method        = PLM + minmod limiter
');
fprintf(fid, 'Time integrator       = RK2
');
fprintf(fid, 'Positivity floors     = enabled
');
fprintf(fid, 'rho floor             = %.6e
', positivity.rhoFloor);
fprintf(fid, 'p floor               = %.6e
', positivity.pFloor);
fprintf(fid, 'Nx = %d
', Nx);
fprintf(fid, 'Ny = %d
', Ny);
fprintf(fid, 'rho0 = %.6f
', params.rho0);
fprintf(fid, 'p0   = %.6f
', params.p0);
fprintf(fid, 'B0   = %.6f
', params.B0);
fprintf(fid, 'A    = %.6e
', params.A);
fprintf(fid, 'mode = %d
', params.mode);
fprintf(fid, 'k    = %.6f
', params.k);
fprintf(fid, 'vA expected = %.6f
', params.vA);
fprintf(fid, 'period      = %.6f
', params.period);
fprintf(fid, 'final time  = %.6f
', t);
fprintf(fid, 'steps       = %d
', step);
fprintf(fid, 'final L2 error in uz = %.6e
', finalErrUz);
fprintf(fid, 'final L2 error in Bz = %.6e
', finalErrBz);
fprintf(fid, 'measured phase speed = %.6f
', measuredSpeed);
fprintf(fid, 'max |div B| final    = %.6e
', maxAbsDivB);
fprintf(fid, 'rho min final        = %.6e
', stateInfo.rhoMin);
fprintf(fid, 'p min final          = %.6e
', stateInfo.pMin);
fprintf(fid, 'positivity floor events      = %d
', numFloorEvents);
fprintf(fid, 'total density-floor cells    = %d
', numRhoFloorCells);
fprintf(fid, 'total pressure-floor cells   = %d
', numPFloorCells);
fclose(fid);

fprintf('Alfven-wave case finished.
');
fprintf('Expected Alfven speed  : %.6f
', params.vA);
fprintf('Measured phase speed   : %.6f
', measuredSpeed);
fprintf('Final L2 error in u_z  : %.6e
', finalErrUz);
fprintf('Final L2 error in B_z  : %.6e
', finalErrBz);
fprintf('Final max |div B|      : %.6e
', maxAbsDivB);
fprintf('Final rho min          : %.6e
', stateInfo.rhoMin);
fprintf('Final p min            : %.6e
', stateInfo.pMin);

results = struct();
results.grid = grid;
results.params = params;
results.positivity = positivity;
results.tFinal = t;
results.stepCount = step;
results.Ufinal = U;
results.Vfinal = Vfinal;
results.timeHistory = timeHistory;
results.dtHistory = dtHistory;
results.errUzHistory = errUzHistory;
results.errBzHistory = errBzHistory;
results.divBmaxHistory = divBmaxHistory;
results.finalErrUz = finalErrUz;
results.finalErrBz = finalErrBz;
results.expectedSpeed = params.vA;
results.measuredSpeed = measuredSpeed;
results.maxAbsDivB = maxAbsDivB;
results.stateInfo = stateInfo;
results.numFloorEvents = numFloorEvents;
results.numRhoFloorCells = numRhoFloorCells;
results.numPFloorCells = numPFloorCells;

end
