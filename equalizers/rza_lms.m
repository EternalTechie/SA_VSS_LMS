function result = rza_lms(data, cfg, alg)
%RZA_LMS Reweighted Zero-Attracting LMS adaptive FIR channel identifier.

N = cfg.N;
M = cfg.M;

mu    = alg.mu;
rho   = alg.rho;
alpha = alg.alpha;

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

    % RZA-LMS zero-attraction term
    attraction = sign(w) ./ (1 + alpha * abs(w));

    % RZA-LMS coefficient update
    w = w + mu * e(n) * dvec - rho * attraction;

    % Store coefficient trajectory
    W(:,n) = w;
end

result.y = y;
result.e = e;
result.w = W;
result.name = "rza_lms";

end