function test_positivity_floor_application()
% test_positivity_floor_application
%
% Verifies that density and pressure floors are applied correctly by the
% positivity-enforcement helper.

clc;

gamma = 5/3;
positivity = default_positivity_settings();
positivity.rhoFloor = 1.0e-6;
positivity.pFloor   = 2.0e-6;

Ny = 2;
Nx = 3;
nVar = 8;

V = zeros(Ny, Nx, nVar);
V(:,:,1) = 1.0;
V(:,:,2) = 0.1;
V(:,:,3) = 0.0;
V(:,:,4) = 0.0;
V(:,:,5) = 1.0;
V(:,:,6) = 0.2;
V(:,:,7) = 0.0;
V(:,:,8) = 0.0;

U = primitive_to_conserved(V, gamma);

% Force one negative density and one negative pressure state.
Vbad = conserved_to_primitive(U, gamma);
Vbad(1,1,1) = -1.0e-3;
Vbad(2,3,5) = -5.0e-4;
Ubad = primitive_to_conserved(Vbad, gamma);

[Ufixed, floorInfo] = apply_positivity_floors(Ubad, gamma, positivity);
Vfixed = conserved_to_primitive(Ufixed, gamma);

rhoMin = min(Vfixed(:,:,1), [], 'all');
pMin   = min(Vfixed(:,:,5), [], 'all');

assert(rhoMin >= positivity.rhoFloor, 'FAIL: Density floor was not enforced.');
assert(pMin   >= positivity.pFloor,   'FAIL: Pressure floor was not enforced.');
assert(floorInfo.rhoFlooredCount >= 1, 'FAIL: Expected at least one density-floor event.');
assert(floorInfo.pFlooredCount   >= 1, 'FAIL: Expected at least one pressure-floor event.');

fprintf('PASS: Positivity floors are applied correctly.\n');

end
