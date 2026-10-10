
function result = mvc(data, cfg, alg)
% MVC - Maximum Versoria Criterion adaptive filter
%
% Required inputs:
%   data.d : Reference input, N x 1
%   data.x : Desired channel output, N x 1
%   cfg.N  : Number of samples
%   cfg.M  : Adaptive filter length
%   alg.mu : Step size
%
% Optional:
%   alg.tau : Versoria parameter (default 1)
%   alg.p   : Versoria exponent (default 2)

    N = cfg.N;
    M = cfg.M;

    mu = alg.mu;

    if isfield(alg, 'tau')
        tau = alg.tau;
    else
        tau = 1;
    end

    if isfield(alg, 'p')
        p = alg.p;
    else
        p = 2;
    end

    if tau <= 0 || p < 1
        error('MVC requires tau > 0 and p >= 1.');
    end

    w = zeros(M, 1);
    W = zeros(M, N);
    y = zeros(N, 1);
    e = zeros(N, 1);

    for n = M:N
        xvec = data.d(n:-1:n-M+1);

        y(n) = w' * xvec;
        e(n) = data.x(n) - y(n);

        % Maximum Versoria score
        score = p * tau * abs(e(n))^(p-1) ...
            * sign(e(n)) ...
            / (1 + tau * abs(e(n))^p)^2;

        w = w + mu * score * xvec;

        W(:, n) = w;
    end

    result.y = y;
    result.e = e;
    result.w = W;
    result.name = 'mvc';
end
