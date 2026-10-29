function result = gsrs(data, cfg, alg)
% GSRS - Generalized Soft-Root-Sign adaptive filter
% (Patel, Bhattacharjee, Christensen, IEEE SPL 2023, doi 10.1109/LSP.2023.3252412)
%
% Cost, with u = |e|^gamma:
%   J(e) = u / (u/alpha + exp(-beta*u))
% Score (eq. 4 of the paper), written in an overflow-safe form:
%   xi = exp(-beta*u) * |e|^(gamma-1) * sign(e) * (beta*u + 1) / (u + alpha*exp(-beta*u))^2
% Update:
%   w = w + mu * xi * x
%
% For small errors and gamma = 2, xi ~ e / alpha^2, so the effective step
% is mu / alpha^2. For large errors xi decays to zero (outlier rejection).
%
% Required inputs:
%   data.d : Reference input, N x 1
%   data.x : Desired channel output, N x 1
%   cfg.N  : Number of samples
%   cfg.M  : Adaptive filter length
%   alg.mu : Step size
%
% Optional:
%   alg.alpha : Scaling parameter (default 3)
%   alg.beta  : Slope parameter (default 1e-3)
%   alg.gamma : Shape parameter (default 2)

    N = cfg.N;
    M = cfg.M;

    mu = alg.mu;

    if isfield(alg, 'alpha')
        alpha = alg.alpha;
    else
        alpha = 3;
    end

    if isfield(alg, 'beta')
        beta = alg.beta;
    else
        beta = 1e-3;
    end

    if isfield(alg, 'gamma')
        gamma = alg.gamma;
    else
        gamma = 2;
    end

    if alpha <= 0 || beta < 0 || gamma <= 0
        error('GSRS requires alpha > 0, beta >= 0 and gamma > 0.');
    end

    w = zeros(M, 1);
    W = zeros(M, N);
    y = zeros(N, 1);
    e = zeros(N, 1);

    for n = M:N
        xvec = data.d(n:-1:n-M+1);

        y(n) = w' * xvec;
        e(n) = data.x(n) - y(n);

        % Generalized soft-root-sign score
        ae = abs(e(n));
        if ae == 0
            score = 0;
        else
            u  = ae^gamma;
            ex = exp(-beta * u);
            score = ex * ae^(gamma - 1) * sign(e(n)) * (beta * u + 1) ...
                    / (u + alpha * ex)^2;
        end

        w = w + mu * score * xvec;

        W(:, n) = w;
    end

    result.y = y;
    result.e = e;
    result.w = W;
    result.name = 'gsrs';
end