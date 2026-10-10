function result = ipnlms(data, cfg, alg)
%IPNLMS Improved Proportionate NLMS adaptive FIR channel identifier
% (Benesty, Gay 2002).
%
% Gains blend NLMS (uniform) and PNLMS (proportionate):
%   k_l = (1 - alpha)/(2M) + (1 + alpha)*|w_l| / (2*||w||_1 + epsK)
%   w   = w + mu * e * K*x / (x'*K*x + delta)
% alpha = -1 gives NLMS; alpha close to 1 behaves like PNLMS.

N = cfg.N;
M = cfg.M;

mu    = alg.mu;
delta = alg.delta;
alpha = alg.alpha;
epsK  = alg.epsK;

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

    % Proportionate gains
    k = (1 - alpha) / (2 * M) + ...
        (1 + alpha) * abs(w) / (2 * sum(abs(w)) + epsK);

    % IPNLMS coefficient update
    kx = k .* dvec;
    w = w + mu * e(n) * kx / (dvec' * kx + delta);

    % Store coefficient trajectory
    W(:, n) = w;

end

result.y = y;
result.e = e;
result.w = W;
result.name = "ipnlms";

end
