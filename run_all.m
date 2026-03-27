function run_all()
% run_all
%
% Convenience script for the current project stage.

addpath(genpath(pwd));

test_variable_conversion;
test_flux_functions;
test_time_step_and_wave_speeds;
test_uniform_state_preservation;
test_rk2_uniform_state_preservation;
test_plm_uniform_state_preservation;

alfven_wave_case;
alfven_wave_convergence_study;
magnetic_pressure_pulse_case;

end
