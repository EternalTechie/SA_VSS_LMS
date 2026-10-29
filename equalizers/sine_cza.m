
function result = sine_cza(data, cfg, alg)
% SINE_CZA - Sine-error LMS with cosine zero attraction (experimental)
%
% Error function (no clamp, by design for this test):
%   f(e) = (T/pi) * sin(pi*e/T)
%   - slope 1 at e = 0, so small errors give the LMS update mu*e*x
%   - rises to T/pi at |e| = T/2, falls to 0 at |e| = T
%   - NEGATIVE for T < |e| < 2T (update points the wrong way), then repeats
%   - always bounded: |f(e)| <= T/pi
%
% Zero attraction (CZA):
%   c(w) = sign(w) * 0.5 * (1 + cos(pi*|w|/a)),  |w| <= a
%          0,                                    |w| >  a
%
% Update:
%   w = w + mu * f(e) * x - rho * c(w)
%   rho = 0 tests the sine error function on its own.
%
% Required inputs:
%   data.d : Reference input, N x 1
%   data.x : Desired channel output, N x 1
%   cfg.N  : Number of samples
%   cfg.M  : Adaptive filter length
%   alg.mu : Step size
%   alg.T  : Error scale; f(e) first returns to zero at |e| = T
%
% Optional:
%   alg.rho : Zero-attraction strength (default 0)
%   alg.a   : Zero-attraction cutoff   (default 0.1)

    N = cfg.N;
    M = cfg.M;

    mu = alg.mu;
    T  = alg.T;

    if isfield(alg, 'rho')
        rho = alg.rho;
    else
        rho = 0;
    end

    if isfield(alg, 'a')
        a = alg.a;
    else
        a = 0.1;
    end

    if T <= 0 || a <= 0 || rho < 0
        error('sine_cza requires T > 0, a > 0 and rho >= 0.');
    end

    w = zeros(M, 1);
    W = zeros(M, N);
    y = zeros(N, 1);
    e = zeros(N, 1);

    for n = M:N
        xvec = data.d(n:-1:n-M+1);

        y(n) = w' * xvec;
        e(n) = data.x(n) - y(n);

        % Sine error function, no clamp
        f = (T / pi) * sin(pi * e(n) / T);

        % Cosine zero attraction
        c = zeros(M, 1);
        idx = abs(w) <= a;
        c(idx) = sign(w(idx)) .* 0.5 .* (1 + cos(pi * abs(w(idx)) / a));

        w = w + mu * f * xvec - rho * c;

        W(:, n) = w;
    end

    result.y = y;
    result.e = e;
    result.w = W;
    result.name = 'sine_cza';
end
