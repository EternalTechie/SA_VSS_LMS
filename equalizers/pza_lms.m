function result = pza_lms(data, cfg, alg)
%PZA_LMS Polynomial Zero-Attracting LMS adaptive FIR channel identifier.

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

    % Polynomial zero-attraction term
    g = zeros(M,1);

    for i = 1:M

        if abs(w(i)) <= 1/(alpha - 1)

            g(i) = sign(w(i)) * ...
                (1 - (alpha - 1)*abs(w(i))) / ...
                (1 + abs(w(i)))^(alpha + 1);

        end

    end

    % PZA-LMS coefficient update
    w = w + mu * e(n) * dvec - rho * g;

    % Store coefficient trajectory
    W(:,n) = w;

end

result.y = y;
result.e = e;
result.w = W;
result.name = "pza_lms";

end