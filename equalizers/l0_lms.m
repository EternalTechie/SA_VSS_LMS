function result = l0_lms(data, cfg, alg)
%L0_LMS l0-norm constrained LMS adaptive FIR channel identifier
% (Gu, Jin, Mei 2009).
%
% Attraction (scaled to peak at 1, like PZA and CZA):
%   c(w) = sign(w) * (1 - beta*|w|),   |w| <= 1/beta
%          0,                          |w| >  1/beta
% rho here equals kappa*beta in the paper's notation.

N = cfg.N;
M = cfg.M;

mu   = alg.mu;
rho  = alg.rho;
beta = alg.beta;

w = zeros(M, 1);

y = zeros(N, 1);
e = zeros(N, 1);
W = zeros(M, N);

for n = M:N

    % Current input vector, newest sample first
    dvec = data.d(n:-1:n-M+1);

    % Estimated channel output
    y(n) = w' * dvec;

    % Desired signal is the observed channel output
    x_des = data.x(n);

    % Instantaneous error
    e(n) = x_des - y(n);

    % l0 zero-attraction term (linear taper)
    c = zeros(M, 1);
    idx = abs(w) <= 1/beta;
    c(idx) = sign(w(idx)) .* (1 - beta * abs(w(idx)));

    % l0-LMS coefficient update
    w = w + mu * e(n) * dvec - rho * c;

    % Store coefficient trajectory
    W(:, n) = w;

end

result.y = y;
result.e = e;
result.w = W;
result.name = "l0_lms";

end
