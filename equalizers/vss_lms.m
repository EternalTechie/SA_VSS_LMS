function result = vss_lms(data, cfg, alg)
%VSS_LMS Variable Step-Size Least Mean Squares channel identifier.
%
% Step-size rule:
%   mu(n+1) = alpha*mu(n) + gamma*e(n)^2
% with lower and upper bounds. This is a baseline variable-step-size method;
% it does not include sparsity-aware regularization.

N = cfg.N;
M = cfg.M;

mu = getAlgField(alg, "mu0", getAlgField(alg, "mu", 0.01));
alpha = getAlgField(alg, "alpha", 0.97);
gamma = getAlgField(alg, "gamma", 1e-3);
muMin = getAlgField(alg, "muMin", 1e-5);
muMax = getAlgField(alg, "muMax", 0.1);

w = zeros(M, 1);

y = zeros(N, 1);
e = zeros(N, 1);
W = zeros(M, N);
muTrace = nan(N, 1);

for n = M:N

    % Current input vector, newest sample first
    dvec = data.d(n:-1:n-M+1);

    % Estimated channel output
    y(n) = w' * dvec;

    % Desired signal is the observed channel output
    x_des = data.x(n);

    % Instantaneous error
    e(n) = x_des - y(n);

    % LMS coefficient update with the current variable step size
    w = w + mu * e(n) * dvec;

    % Store coefficient and step-size trajectories
    W(:, n) = w;
    muTrace(n) = mu;

    % Error-power controlled variable step-size update
    mu = alpha * mu + gamma * e(n)^2;
    mu = min(max(mu, muMin), muMax);
end

result.y = y;
result.e = e;
result.w = W;
result.mu = muTrace;
result.name = "vss_lms";
end

function value = getAlgField(alg, fieldName, defaultValue)
if isfield(alg, fieldName)
    value = alg.(fieldName);
else
    value = defaultValue;
end
end
