%% Impulse Inspection
% Shows what the impulsive noise looks like and how each robust error
% function treats errors of that size.
% Calls main.m once to initialize cfg and alg, then rebuilds the same
% received signal main.m uses (same seed).

clear; clc; close all;

%% Initialize the main simulation
setappdata(0, 'SAVSS_INIT_ONLY', true);
run('main.m');
if isappdata(0, 'SAVSS_INIT_ONLY')
    rmappdata(0, 'SAVSS_INIT_ONLY');
end
close all;

%% Settings
insp.window     = [1 5000];      % sample range for the time plots
insp.errorAlg   = "andrew_sine"; % filter whose error signal is shown
insp.eZoom      = 5;             % error range for the zoomed error-function plot
insp.eWide      = 80;            % error range for the wide error-function plot

%% Rebuild the received signal exactly as main.m does
rng(cfg.seed);
d    = cfg.signal.generator(cfg.N, cfg.signal);
data = generateReceivedData(d, cfg);

xClean = applyTimeVaryingFIR(d, data.h);   % noiseless channel output
noise  = data.x(:) - xClean(:);            % AWGN + impulses
xClean = xClean(:);

sigPow   = mean(xClean.^2);
sigmaV   = sqrt(sigPow / 10^(cfg.noise.SNR_dB / 10));   % background noise std
isImp    = abs(noise) > 6 * sigmaV;                       % samples dominated by an impulse
impMag   = abs(noise(isImp));

%% Print impulse statistics
fprintf("\n=== Received signal ===\n");
fprintf("Signal RMS            : %.3f\n", sqrt(sigPow));
fprintf("Background noise std  : %.4f  (%g dB SNR)\n", sigmaV, cfg.noise.SNR_dB);

if isfield(cfg.noise, "factors") && ~isempty(cfg.noise.factors)
    fprintf("Noise factors         : %s\n", ...
            strjoin(cellfun(@func2str, cfg.noise.factors, 'UniformOutput', false), ", "));
else
    fprintf("Noise factors         : none\n");
end

fprintf("\n=== Impulses (|noise| > 6 x background std) ===\n");
fprintf("Count                 : %d of %d samples (%.2f%%)\n", ...
        nnz(isImp), numel(noise), 100 * mean(isImp));
if ~isempty(impMag)
    q = prctile(impMag, [10 50 90 99]);
    fprintf("Magnitude 10th / median / 90th / 99th percentile: %.1f / %.1f / %.1f / %.1f\n", q);
    fprintf("Largest impulse       : %.1f  (%.0f x signal RMS, %.0f x background std)\n", ...
            max(impMag), max(impMag) / sqrt(sigPow), max(impMag) / sigmaV);
end

%% Plot styles
grayC   = [0.55 0.55 0.55];
blueC   = [0.16 0.47 0.84];
orangeC = [0.92 0.41 0.20];

w = insp.window(1):min(insp.window(2), cfg.N);

%% Figure 1: received signal and noise over time
figure;
tiledlayout(2, 1);

nexttile;
plot(w, data.x(w), 'Color', orangeC);
hold on;
plot(w, xClean(w), 'Color', blueC, 'LineWidth', 1);
hold off;
grid on;
xlabel('Sample');
ylabel('Amplitude');
title('Received signal x (orange) vs clean channel output (blue)');

nexttile;
plot(w, noise(w), 'Color', grayC);
hold on;
idx = w(isImp(w));
plot(idx, noise(idx), 'o', 'Color', orangeC, 'MarkerSize', 5);
hold off;
grid on;
xlabel('Sample');
ylabel('Noise');
title('Noise = AWGN + impulses (circles = impulses)');

%% Figure 2: impulse magnitude distribution
if ~isempty(impMag)
    figure;
    histogram(impMag, 40, 'FaceColor', orangeC);
    grid on;
    xlabel('Impulse magnitude |v|');
    ylabel('Count');
    title(sprintf('Impulse magnitudes (%d impulses, median %.1f)', ...
          numel(impMag), median(impMag)));
end

%% Run one filter to get a real error signal
res = runEqualizer(insp.errorAlg, data, cfg, alg);
e   = res.e(:);
eSS = e(floor(0.5 * cfg.N):end);           % second half, after convergence

%% Error functions, each scaled to the common step size
% g(e) = mu_alg * f(e) / mu_common, so every curve has slope 1 at e = 0
funcs = buildErrorFunctions(alg, mu_common);

for view = 1:2
    if view == 1
        eMax = insp.eZoom;  yl = [-1.5 1.5] * insp.eZoom / 2;
        ttl = 'Error functions near zero (normal errors)';
    else
        eMax = insp.eWide;  yl = [-12 12];
        ttl = 'Error functions over the impulse range';
    end
    eg = linspace(-eMax, eMax, 4001);

    figure;
    yyaxis right;
    histogram(eSS(abs(eSS) <= eMax), linspace(-eMax, eMax, 121), ...
              'FaceColor', grayC, 'EdgeColor', 'none', 'FaceAlpha', 0.35, ...
              'DisplayName', "errors e(n)");
    set(gca, 'YScale', 'log');
    ylabel(sprintf('Count of errors e(n) from %s (log scale)', insp.errorAlg), ...
           'Interpreter', 'none');

    yyaxis left;
    hold on;
    colors     = turbo(numel(funcs) + 2);
    colors     = colors(2:end-1, :);
    lineStyles = ["-", "--", "-.", ":"];
    for k = 1:numel(funcs)
        plot(eg, funcs(k).g(eg), ...
            'Color', colors(k, :), ...
            'LineStyle', lineStyles(mod(k - 1, numel(lineStyles)) + 1), ...
            'LineWidth', 1.6, 'Marker', 'none', ...
            'DisplayName', funcs(k).name);
    end
    hold off;
    ylim(yl);
    ylabel('Update strength g(e) (LMS = e)');
    xlabel('Error e');
    grid on;
    title(ttl);
    legend('Interpreter', 'none', 'Location', 'eastoutside');
end

%% Figure: error signal over time
figure;
plot(e, 'Color', grayC);
grid on;
xlabel('Sample');
ylabel('e(n)');
title(sprintf('Error signal of %s (spikes are impulses)', insp.errorAlg), ...
      'Interpreter', 'none');

%% Local functions

function funcs = buildErrorFunctions(alg, muC)
%BUILDERRORFUNCTIONS Error functions of each filter, scaled to muC.

funcs = struct('name', {}, 'g', {});

funcs(end+1) = struct('name', "lms", 'g', @(e) e);

if isfield(alg, 'mvc')
    t = alg.mvc.tau;  m = alg.mvc.mu / muC;
    funcs(end+1) = struct('name', "mvc", 'g', @(e) m * 2 * t * e ./ (1 + t * e.^2).^2);
end
if isfield(alg, 'andrew_sine')
    c = alg.andrew_sine.c;  m = alg.andrew_sine.mu / muC;
    funcs(end+1) = struct('name', "andrew_sine", ...
                          'g', @(e) m * sin(e / c) .* (abs(e) <= pi * c));
end
if isfield(alg, 'hcaf')
    l = alg.hcaf.lambda;  m = alg.hcaf.mu / muC;
    funcs(end+1) = struct('name', "hcaf", 'g', @(e) m * sinh(l * e));
end
if isfield(alg, 'lhcaf')
    l = alg.lhcaf.lambda;  m = alg.lhcaf.mu / muC;
    funcs(end+1) = struct('name', "lhcaf / jlhcaf", 'g', @(e) m * tanh(l * e));
end
if isfield(alg, 'gsrs')
    a = alg.gsrs.alpha;  b = alg.gsrs.beta;  m = alg.gsrs.mu / muC;
    funcs(end+1) = struct('name', "gsrs", 'g', @(e) m * gsrsScore(e, a, b));
end
if isfield(alg, 'sine_cza')
    T = alg.sine_cza.T;  m = alg.sine_cza.mu / muC;
    funcs(end+1) = struct('name', "sine_cza", 'g', @(e) m * (T / pi) * sin(pi * e / T));
end
end

function s = gsrsScore(e, alpha, beta)
%GSRSSCORE GSRS score for gamma = 2 (overflow-safe form).
u  = e.^2;
ex = exp(-beta * u);
s  = ex .* e .* (beta * u + 1) ./ (u + alpha * ex).^2;
end