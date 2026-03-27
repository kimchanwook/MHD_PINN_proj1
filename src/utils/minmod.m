function mm = minmod(a, b)
% minmod
%
% Elementwise minmod limiter for piecewise-linear reconstruction.
%
% INPUTS:
%   a - first slope estimate
%   b - second slope estimate
%
% OUTPUT:
%   mm - minmod-limited slope
%
% DEFINITION:
%   If a and b have the same sign, the result is the one with smaller
%   magnitude. If they have opposite sign, the result is zero.

mm = zeros(size(a));
sameSign = (a .* b) > 0;
mm(sameSign) = sign(a(sameSign)) .* min(abs(a(sameSign)), abs(b(sameSign)));

end
