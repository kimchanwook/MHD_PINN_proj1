function errL2 = compute_l2_error(numericalField, exactField)
% compute_l2_error
%
% Computes the discrete L2 error norm between a numerical field and an exact
% reference field on the same grid.

diffField = numericalField - exactField;
errL2 = sqrt(mean(diffField(:).^2));

end
