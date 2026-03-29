function results = magnetic_pressure_pulse_case()
% magnetic_pressure_pulse_case
%
% Runs a nonlinear magnetic-pressure pulse demonstration using second-order
% Runge-Kutta time stepping and piecewise-linear reconstruction.
%
% This stage also enables optional positivity floors for density and gas
% pressure to improve robustness in strong-gradient regions.

clc;

gamma = 5/3;
CFL   = 0.35;
positivity = default_positivity_settings();

Nx = 128;
Ny = 128;

xMin = 0.0;
xMax = 1.0;
yMin = 0.0;
yMax = 1.0;

gridData = make_uniform_grid(xMin, xMax, yMin, yMax, Nx, Ny);

caseParams = struct();
caseParams.rho0  = 1.0;
caseParams.p0    = 1.0;
caseParams.Bx0   = 1.0;
caseParams.By0   = 0.0;
caseParams.A     = 0.5;
caseParams.sigma = 0.08;
caseParams.x0    = 0.5;
caseParams.y0    = 0.5;

outDir = fullfile('output', 'magnetic_pressure_pulse');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

[U, params] = init_magnetic_pressure_pulse(gridData, gamma, caseParams);

t = 0.0;
finalTime = 0.2;
step = 0;
maxSteps = 300000;

sampleTimes = [0.0, 0.05, 0.10, 0.20];
sampleIndex = 1;
numFloorEvents = 0;
numRhoFloorCells = 0;
numPFloorCells = 0;

massHistory = [];
energyHistory = [];
divBmaxHistory = [];
timeHistory = [];

V = conserved_to_primitive(U, gamma);
plot_pulse_fields(gridData, V, t, outDir, 'pulse_t0');

while t < finalTime && step < maxSteps
    dt = compute_time_step(U, gridData, gamma, CFL);
    if t + dt > finalTime
        dt = finalTime - t;
    end

    [U, stepInfo] = update_fv_rk2_plm(U, gridData, gamma, dt, positivity);
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

    divB = compute_divB(U, gridData);
    massHistory(end+1,1) = compute_total_mass(U, gridData); %#ok<AGROW>
    energyHistory(end+1,1) = compute_total_energy_domain(U, gridData); %#ok<AGROW>
    divBmaxHistory(end+1,1) = max(abs(divB(:))); %#ok<AGROW>
    timeHistory(end+1,1) = t; %#ok<AGROW>

    while sampleIndex <= numel(sampleTimes) && t >= sampleTimes(sampleIndex) - 1e-12
        Vsample = conserved_to_primitive(U, gamma);
        prefix = sprintf('pulse_t%03d', round(1000*sampleTimes(sampleIndex)));
        plot_pulse_fields(gridData, Vsample, t, outDir, prefix);
        sampleIndex = sampleIndex + 1;
    end
end

Vfinal = conserved_to_primitive(U, gamma);
stateInfo = check_physical_state(U, gamma);
plot_pulse_fields(gridData, Vfinal, t, outDir, 'pulse_final');

fig1 = figure('Visible','off');
plot(timeHistory, massHistory, 'LineWidth', 1.5);
xlabel('time'); ylabel('Total mass');
title('Magnetic-pressure pulse: total mass');
grid on;
saveas(fig1, fullfile(outDir, 'pulse_mass_history.png'));
close(fig1);

fig2 = figure('Visible','off');
plot(timeHistory, energyHistory, 'LineWidth', 1.5);
xlabel('time'); ylabel('Total energy');
title('Magnetic-pressure pulse: total energy');
grid on;
saveas(fig2, fullfile(outDir, 'pulse_energy_history.png'));
close(fig2);

fig3 = figure('Visible','off');
plot(timeHistory, divBmaxHistory, 'LineWidth', 1.5);
xlabel('time'); ylabel('max |div B|');
title('Magnetic-pressure pulse: divergence diagnostic');
grid on;
saveas(fig3, fullfile(outDir, 'pulse_divB_history.png'));
close(fig3);

fid = fopen(fullfile(outDir, 'pulse_summary.txt'), 'w');
fprintf(fid, 'Magnetic pressure pulse summary
');
fprintf(fid, '--------------------------------
');
fprintf(fid, 'Spatial method  = PLM + minmod limiter
');
fprintf(fid, 'Time integrator = RK2
');
fprintf(fid, 'Positivity floors = enabled
');
fprintf(fid, 'rho floor = %.6e
', positivity.rhoFloor);
fprintf(fid, 'p floor   = %.6e
', positivity.pFloor);
fprintf(fid, 'Nx = %d
', Nx);
fprintf(fid, 'Ny = %d
', Ny);
fprintf(fid, 'rho0  = %.6f
', params.rho0);
fprintf(fid, 'p0    = %.6f
', params.p0);
fprintf(fid, 'Bx0   = %.6f
', params.Bx0);
fprintf(fid, 'By0   = %.6f
', params.By0);
fprintf(fid, 'A     = %.6f
', params.A);
fprintf(fid, 'sigma = %.6f
', params.sigma);
fprintf(fid, 'x0    = %.6f
', params.x0);
fprintf(fid, 'y0    = %.6f
', params.y0);
fprintf(fid, 'final time = %.6f
', t);
fprintf(fid, 'steps      = %d
', step);
fprintf(fid, 'final total mass   = %.12e
', massHistory(end));
fprintf(fid, 'final total energy = %.12e
', energyHistory(end));
fprintf(fid, 'final max |div B|  = %.6e
', divBmaxHistory(end));
fprintf(fid, 'final rho min      = %.6e
', stateInfo.rhoMin);
fprintf(fid, 'final p min        = %.6e
', stateInfo.pMin);
fprintf(fid, 'positivity floor events      = %d
', numFloorEvents);
fprintf(fid, 'total density-floor cells    = %d
', numRhoFloorCells);
fprintf(fid, 'total pressure-floor cells   = %d
', numPFloorCells);
fclose(fid);

results = struct();
results.gridData = gridData;
results.params = params;
results.positivity = positivity;
results.tFinal = t;
results.stepCount = step;
results.Ufinal = U;
results.Vfinal = Vfinal;
results.timeHistory = timeHistory;
results.massHistory = massHistory;
results.energyHistory = energyHistory;
results.divBmaxHistory = divBmaxHistory;
results.stateInfo = stateInfo;
results.numFloorEvents = numFloorEvents;
results.numRhoFloorCells = numRhoFloorCells;
results.numPFloorCells = numPFloorCells;

fprintf('Magnetic-pressure pulse case finished.
');
fprintf('Final total mass   : %.12e
', massHistory(end));
fprintf('Final total energy : %.12e
', energyHistory(end));
fprintf('Final max |div B|  : %.6e
', divBmaxHistory(end));
fprintf('Final rho min      : %.6e
', stateInfo.rhoMin);
fprintf('Final p min        : %.6e
', stateInfo.pMin);

end
