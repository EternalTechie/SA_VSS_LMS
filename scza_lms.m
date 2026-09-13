function result = scza_lms(data, cfg, alg)
%SATURATING_COSINE_LMS LMS adaptive FIR channel identifier with
% saturating-cosine zero attraction.

N = cfg.N;
M = cfg.M;

mu = alg.mu;
rho = alg.rho;
a  = alg.a;

w = zeros(M,1);

y = zeros(N,1);
e = zeros(N,1);
W = zeros(M,N);

for n = M:N

    % Current input vector, newest sample first
    dvec = data.d(n:-1:n-M+1);

    % Estimated channel output
    y(n) = w' * dvec;

    % Desired signal is the observed channel output
    x_des = data.x(n);

    % Instantaneous error
    e(n) = x_des - y(n);

    % Saturating-cosine zero-attraction
    attraction = sign(w) .* cos( ...
        (pi/2) .* abs(w) ./ (a + abs(w)) );

    % Coefficient update
    w = w + mu * e(n) * dvec - rho * attraction;

    % Store coefficient trajectory
    W(:,n) = w;

end

result.y = y;
result.e = e;
result.w = W;
result.name = "scza_lms";

end