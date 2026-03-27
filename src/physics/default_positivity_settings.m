function positivity = default_positivity_settings()
% default_positivity_settings
%
% Returns a default positivity-floor configuration for the ideal-MHD solver.

positivity = struct();
positivity.isEnabled = true;
positivity.rhoFloor  = 1.0e-8;
positivity.pFloor    = 1.0e-8;

end
