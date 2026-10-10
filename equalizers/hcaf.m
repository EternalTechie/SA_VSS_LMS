
function result = hcaf(data, cfg, alg)
% HCAF - Hyperbolic Cosine Adaptive Filter
%
% Required:
%   alg.mu     : Step size
%
% Optional:
%   alg.lambda : Hyperbolic scaling parameter (default 1)

    N = cfg.N;
    M = cfg.M;

    mu = alg.mu;

    if isfield(alg, 'lambda')
        lambda = alg.lambda;
    else
        lambda = 1;
    end

    if lambda <= 0
        error('HCAF requires lambda > 0.');
    end

    w = zeros(M, 1);
    W = zeros(M, N);
    y = zeros(N, 1);
    e = zeros(N, 1);

    for n = M:N
        xvec = data.d(n:-1:n-M+1);

        y(n) = w' * xvec;
        e(n) = data.x(n) - y(n);

        % Hyperbolic cosine score
        score = sinh(lambda * e(n));

        w = w + mu * score * xvec;

        W(:, n) = w;
    end

    result.y = y;
    result.e = e;
    result.w = W;
    result.name = 'hcaf';
end
