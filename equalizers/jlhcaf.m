
function result = jlhcaf(data, cfg, alg)
% JLHCAF - Joint logarithmic hyperbolic cosine adaptive filter
% (Kumar, Bhattacharjee, George, IEEE TCAS-II, vol. 68, no. 1, 2021, eq. 8)
%
% Cost:
%   xi(n) = ln cosh(lambda*e(n)) + (gamma/a) * sum_j ln cosh(a*w_j(n))
% Update:
%   w = w + mu * tanh(lambda*e) * x - rho * tanh(a*w)
%
% Error term: LMS-like for small errors (slope mu*lambda), capped at mu
% for large errors.
% Weight term: tanh(a*w) is ~a*w near zero (a leakage proportional to w)
% and saturates at 1 for large |a*w| (a constant pull like ZA-LMS).
%
% Required inputs:
%   data.d     : Reference input, N x 1
%   data.x     : Desired channel output, N x 1
%   cfg.N      : Number of samples
%   cfg.M      : Adaptive filter length
%   alg.mu     : Step size
%   alg.lambda : Error scaling factor
%   alg.rho    : Penalty strength
%   alg.a      : Weight scaling factor in tanh(a*w)

    N = cfg.N;
    M = cfg.M;

    mu     = alg.mu;
    lambda = alg.lambda;
    rho    = alg.rho;
    a      = alg.a;

    if lambda <= 0 || a <= 0 || rho < 0
        error('jlhcaf requires lambda > 0, a > 0 and rho >= 0.');
    end

    w = zeros(M, 1);
    W = zeros(M, N);
    y = zeros(N, 1);
    e = zeros(N, 1);

    for n = M:N
        xvec = data.d(n:-1:n-M+1);

        y(n) = w' * xvec;
        e(n) = data.x(n) - y(n);

        w = w + mu * tanh(lambda * e(n)) * xvec - rho * tanh(a * w);

        W(:, n) = w;
    end

    result.y = y;
    result.e = e;
    result.w = W;
    result.name = 'jlhcaf';
end
