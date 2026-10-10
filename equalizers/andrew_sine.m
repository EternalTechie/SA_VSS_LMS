
function result = andrew_sine(data, cfg, alg)
% ANDREW_SINE - Andrew's Sine robust adaptive filter
%
% Required:
%   alg.mu : Step size
%
% Optional:
%   alg.c  : Scale parameter (default 1)

    N = cfg.N;
    M = cfg.M;

    mu = alg.mu;

    if isfield(alg, 'c')
        c = alg.c;
    else
        c = 1;
    end

    if c <= 0
        error('Andrew''s Sine requires c > 0.');
    end

    w = zeros(M, 1);
    W = zeros(M, N);
    y = zeros(N, 1);
    e = zeros(N, 1);

    for n = M:N
        xvec = data.d(n:-1:n-M+1);

        y(n) = w' * xvec;
        e(n) = data.x(n) - y(n);

        % Reject errors outside the specified interval
        if abs(e(n)) <= pi * c
            score = sin(e(n) / c);
            w = w + mu * score * xvec;
        end

        W(:, n) = w;
    end

    result.y = y;
    result.e = e;
    result.w = W;
    result.name = 'andrew_sine';
end
