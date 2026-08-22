function x = applyTimeVaryingFIR(d, h)
%APPLYTIMEVARYINGFIR Apply a possibly time-varying finite impulse response.

N = numel(d);
M = size(h, 1);
x = zeros(N, 1);

for n = 1:N
    tapCount = min(M, n);
    dvec = d(n:-1:n-tapCount+1);
    x(n) = h(1:tapCount, n)' * dvec;
end
end
