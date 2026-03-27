function measuredSpeed = measure_alfven_phase_speed(x, fieldAtT0, fieldAtT, t)
% measure_alfven_phase_speed
%
% Estimates the phase speed of a periodic traveling wave by computing the
% spatial shift that best aligns fieldAtT with fieldAtT0.

x = x(:);
f0 = fieldAtT0(:);
ft = fieldAtT(:);

Nx = length(x);
dx = x(2) - x(1);
Lx = Nx * dx;

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

if shiftDistance > 0.5 * Lx
    shiftDistance = shiftDistance - Lx;
end

measuredSpeed = shiftDistance / t;

end
