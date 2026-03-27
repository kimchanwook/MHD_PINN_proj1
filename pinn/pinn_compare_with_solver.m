function comparison = pinn_compare_with_solver(net, grid, exact, tEval)
% pinn_compare_with_solver
%
% Compares PINN predictions against the analytical Alfvén-wave reference on
% the solver grid at a selected evaluation time.
%
% INPUTS:
%   net   - PINN network
%   grid  - finite-volume grid structure
%   exact - exact Alfvén-wave solution structure
%   tEval - evaluation time
%
% OUTPUT:
%   comparison - structure with predicted and exact line cuts and L2 errors

x = grid.xc(:);
yMid = mean([grid.yMin, grid.yMax]);
N = numel(x);

xyt = [x, yMid * ones(N,1), tEval * ones(N,1)];
pred = pinn_predict_fields(net, xyt);

uzExact = exact.uz(x, tEval);
BzExact = exact.Bz(x, tEval);

comparison = struct();
comparison.x = x;
comparison.tEval = tEval;
comparison.rhoPred = pred.rho;
comparison.uxPred = pred.ux;
comparison.uyPred = pred.uy;
comparison.uzPred = pred.uz;
comparison.pPred = pred.p;
comparison.BxPred = pred.Bx;
comparison.ByPred = pred.By;
comparison.BzPred = pred.Bz;
comparison.uzExact = uzExact;
comparison.BzExact = BzExact;
comparison.errUzL2 = compute_l2_error(pred.uz, uzExact);
comparison.errBzL2 = compute_l2_error(pred.Bz, BzExact);
comparison.maxAbsRhoDrift = max(abs(pred.rho - exact.rho0));
comparison.maxAbsPDrift   = max(abs(pred.p - exact.p0));
comparison.maxAbsBxDrift  = max(abs(pred.Bx - exact.B0));
comparison.maxAbsUy       = max(abs(pred.uy));
comparison.maxAbsBy       = max(abs(pred.By));

end
