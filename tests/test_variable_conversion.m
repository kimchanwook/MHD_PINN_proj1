function test_variable_conversion()
% test_variable_conversion
%
% Basic consistency test for primitive/conserved conversion.

clc;

gamma = 5/3;

Ny = 4;
Nx = 5;
nVar = 8;

V = zeros(Ny, Nx, nVar);

V(:,:,1) = 1.2;
V(:,:,2) = 0.3;
V(:,:,3) = -0.1;
V(:,:,4) = 0.05;
V(:,:,5) = 0.9;
V(:,:,6) = 0.7;
V(:,:,7) = -0.2;
V(:,:,8) = 0.15;

U = primitive_to_conserved(V, gamma);
V_recovered = conserved_to_primitive(U, gamma);

abs_error = abs(V_recovered - V);
max_abs_error = max(abs_error(:));

fprintf('Max absolute error in V -> U -> V recovery: %.16e\n', max_abs_error);

tolerance = 1e-12;

if max_abs_error < tolerance
    fprintf('PASS: Variable conversion test passed.\n');
else
    error('FAIL: Variable conversion test failed.');
end

end
