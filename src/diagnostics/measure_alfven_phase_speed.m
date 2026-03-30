function measuredSpeed = measure_alfven_phase_speed(x, fieldAtT0, fieldAtT, t)
% measure_alfven_phase_speed
%
% Estimates the phase speed of a periodic traveling wave by computing the
% circular shift that best aligns fieldAtT with fieldAtT0.

x = x(:);
f0 = fieldAtT0(:);
ft = fieldAtT(:);

if numel(x) < 2
    error('x must contain at least two points.');
end

if numel(f0) ~= numel(x) || numel(ft) ~= numel(x)
    error('x, fieldAtT0, and fieldAtT must have the same length.');
end

Nx = length(x);
dx = x(2) - x(1);
Lx = Nx * dx;

% Remove any tiny mean offset before correlation
f0 = f0 - mean(f0);
ft = ft - mean(ft);

bestScore = -inf;
bestShiftIndex = 0;

for shiftIndex = 0:(Nx - 1)
    shifted = circshift(f0, shiftIndex);
    score = sum(shifted .* ft);
    if score > bestScore
        bestScore = score;
        bestShiftIndex = shiftIndex;
    end
end

shiftDistance = bestShiftIndex * dx;

% Map shift into [-Lx/2, Lx/2]
if shiftDistance > 0.5 * Lx
    shiftDistance = shiftDistance - Lx;
end

measuredSpeed = shiftDistance / t;

end