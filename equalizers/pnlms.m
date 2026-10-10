function result = pnlms(data, cfg, alg)
%PNLMS Proportionate NLMS adaptive FIR channel identifier (Duttweiler 2000).
%
% Each tap gets a step size proportional to its current magnitude:
%   gamma_l = max(rhoP * max(deltaP, ||w||_inf), |w_l|)
%   g_l     = gamma_l / mean(gamma)
%   w       = w + mu * e * G*x / (x'*G*x + delta)

N = cfg.N;
M = cfg.M;

mu     = alg.mu;
delta  = alg.delta;
rhoP   = alg.rhoP;      % keeps small taps adapting
deltaP = alg.deltaP;    % gets adaptation started from w = 0

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

    % Proportionate gains (average gain = 1)
    gamma = max(rhoP * max(deltaP, max(abs(w))), abs(w));
    g = gamma / mean(gamma);

    % PNLMS coefficient update
    gx = g .* dvec;
    w = w + mu * e(n) * gx / (dvec' * gx + delta);

    % Store coefficient trajectory
    W(:, n) = w;

end

result.y = y;
result.e = e;
result.w = W;
result.name = "pnlms";

end
