function result = cza_lms(data, cfg, alg)

N = cfg.N;
M = cfg.M;

mu  = alg.mu;
rho = alg.rho;
a   = alg.a;

w = zeros(M, 1);

y = zeros(N, 1);
e = zeros(N, 1);
W = zeros(M, N);

for n = M:N

    % Input vector
    x = data.d(n:-1:n-M+1);

    % Equalizer output
    y(n) = w' * x;

    % Desired signal
    d = data.x(n);

    % Error
    e(n) = d - y(n);

    % ---------------------------------------------------------
    % Cosine Zero Attraction LMS
    %
    % c(w) = sign(w)/2 * (1 + cos(pi*|w|/a)), |w| <= a
    %        0,                                  |w| > a
    % ---------------------------------------------------------

    c = zeros(M, 1);

    idx = abs(w) <= a;

    c(idx) = sign(w(idx)) .* ...
             0.5 .* ...
             (1 + cos(pi * abs(w(idx)) / a));

    % LMS + zero-attraction update
    w = w + mu * e(n) * x - rho * c;

    % Store weights
    W(:, n) = w;

end

result.y = y;
result.e = e;
result.w = W;
result.name = "cza_lms";

end