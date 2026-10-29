%% SA-VSS-LMS Adaptive Filtering Simulation
% Main entry point.
% Change the configured signal, channel, noise, and adaptive filter factors
% here to run a different combination under the same evaluation framework.

clear; clc; close all;

addpath(genpath(pwd));

%% Select adaptive filter

equalizerNames = ["rza_lms", ...
                  "pza_lms", "cza_lms", ...
                  "mvc", "andrew_sine", "hcaf", "lhcaf", "sine_cza", "jlhcaf"];

%% Global simulation parameters

cfg.N          = 100000;     % number of samples
cfg.M          = 32;        % adaptive filter length
cfg.SNR_dB     = 30;        % received-signal SNR
cfg.delay      = 0;         % desired-signal alignment

% Discrete-time input/training signal
cfg.signal.generator = @whiteGaussianSignal;

% Channel model
cfg.channel.initializer = @sparseChannel;
cfg.channel.M           = cfg.M;
cfg.channel.activeTaps  = [2 4 9];
cfg.channel.gains       = [0.2 0.8 0.7];
cfg.channel.normalize   = true;
cfg.channel.slowDrift.stdPerSample = 1e-4;
cfg.channel.factors     = {};        % Example: {@slowDriftChannel}

% Receiver noise model
cfg.noise.apply   = @awgnNoise;
cfg.noise.SNR_dB  = cfg.SNR_dB;
cfg.noise.factors = {@impulseNoise};    % Example: {@impulseNoise}

% Impulse settings (used only if @impulseNoise is in cfg.noise.factors)
cfg.noise.impulse.prob       = 0.01;    % 1% of samples hit by impulses
cfg.noise.impulse.powerRatio = 1000;    % impulse variance = 1000 x signal power
cfg.noise.impulse.outlierIdx = [];      % e.g. [20000 40000] for fixed spikes
cfg.noise.impulse.outlierAmp = 100;

% Reproducibility
cfg.seed        = 1;

%% Algorithm parameters
% Common effective step size for fixed-step algorithms
mu_common = 0.01;
M   = 32;      % filter length
sx2 = 1;       % input variance
% Common attraction strength for zero-attracting algorithms
rho_common = 1e-4;
% LMS
alg.lms_c.mu = mu_common;
% NLMS
alg.nlms.mu    = mu_common * M * sx2;   % 0.32, NLMS has input-energy normalization
alg.nlms.delta = 1e-8;
% VSS-LMS
alg.vss_lms.mu0   = mu_common;
alg.vss_lms.alpha = 0.97;
alg.vss_lms.gamma = 1e-3;
alg.vss_lms.muMin = 1e-5;
alg.vss_lms.muMax = mu_common;          % can only adapt down from the common step
% ZA-LMS
alg.za_lms.mu  = mu_common;
alg.za_lms.rho = rho_common;
% RZA-LMS
alg.rza_lms.mu    = mu_common;
alg.rza_lms.rho   = rho_common;
alg.rza_lms.alpha = 10;                 % pull halves at |w| = 0.1
% PZA-LMS
alg.pza_lms.mu    = mu_common;
alg.pza_lms.rho   = rho_common;
alg.pza_lms.alpha = 11;                 % cutoff 1/(alpha-1) = 0.1
% CZA-LMS
alg.cza_lms.mu  = mu_common;
alg.cza_lms.rho = rho_common;
alg.cza_lms.a   = 0.1;                  % cutoff 0.1
% Maximum Versoria Criterion
alg.mvc.tau = 1;
alg.mvc.p   = 2;
alg.mvc.mu  = mu_common / (2 * alg.mvc.tau);         % 0.005, effective step 2*mu*tau
% Andrew's Sine Estimator
alg.andrew_sine.c  = 0.5;
alg.andrew_sine.mu = mu_common * alg.andrew_sine.c;  % 0.005, effective step mu/c
% Log-Hyperbolic Cosine Adaptive Filter
alg.lhcaf.lambda = 2;
alg.lhcaf.mu     = mu_common / alg.lhcaf.lambda;     % 0.005, effective step mu*lambda
% Hyperbolic Cosine Adaptive Filter
alg.hcaf.lambda = 0.1;
alg.hcaf.mu     = mu_common / alg.hcaf.lambda;       % 0.1, effective step mu*lambda
% l0-LMS
alg.l0_lms.mu   = mu_common;
alg.l0_lms.rho  = rho_common;           % rho = kappa*beta in Gu et al. notation
alg.l0_lms.beta = 10;                   % cutoff 1/beta = 0.1
% PNLMS
alg.pnlms.mu     = alg.nlms.mu;         % same average step as NLMS
alg.pnlms.delta  = 1e-8;
alg.pnlms.rhoP   = 0.01;
alg.pnlms.deltaP = 0.01;
% IPNLMS
alg.ipnlms.mu    = alg.nlms.mu;         % same average step as NLMS
alg.ipnlms.delta = 1e-8;
alg.ipnlms.alpha = 0;
alg.ipnlms.epsK  = 1e-6;
% Generalized Soft-Root-Sign (GSRS)
alg.gsrs.alpha = 3;                              % score peaks near |e| = sqrt(alpha/3) = 1
alg.gsrs.beta  = 1e-3;                           % slope parameter
alg.gsrs.gamma = 2;                              % shape; 2 keeps the small-error update linear
alg.gsrs.mu    = mu_common * alg.gsrs.alpha^2;   % 0.09, effective step mu/alpha^2
% Sine-error CZA (experimental, no clamp)
alg.sine_cza.mu  = mu_common;           % slope of f(e) at 0 is 1, so effective step = mu
alg.sine_cza.T   = 3;                   % f(e) returns to zero at |e| = T
alg.sine_cza.rho = rho_common;          % set to 0 to test the sine error function alone
alg.sine_cza.a   = 0.1;                 % same cutoff as CZA-LMS
% Joint Logarithmic Hyperbolic Cosine Adaptive Filter
alg.jlhcaf.lambda = 2;                              % same as LHCAF
alg.jlhcaf.mu     = mu_common / alg.jlhcaf.lambda;  % 0.005, effective step mu*lambda
alg.jlhcaf.rho    = rho_common;
alg.jlhcaf.a      = 0.5;                            % weight scaling in tanh(a*w)

%% Support for ablation_study or monteCarlo etc
if isappdata(0, 'SAVSS_INIT_ONLY') && getappdata(0, 'SAVSS_INIT_ONLY')
    return;
end

%% Generate input signal and received data

rng(cfg.seed);

d = cfg.signal.generator(cfg.N, cfg.signal);
data = generateReceivedData(d, cfg);

%% Run adaptive equalizer

results = struct();

for k = 1:numel(equalizerNames)

    name = equalizerNames(k);

    fprintf("\nRunning %s...\n", name);

    results.(name) = runEqualizer(name, data, cfg, alg);

end

%% Evaluate all equalizers

metrics = struct();

for k = 1:numel(equalizerNames)

    name = equalizerNames(k);

    metrics.(name) = evaluateResult( ...
        data, results.(name), cfg);

end

%% Display basic results

fprintf("\nScenario: %s\n\n", data.name);

for k = 1:numel(equalizerNames)

    name = equalizerNames(k);

    fprintf("%-10s | NMSD: %8.3f dB | MSD: %8.3f dB | MSE: %8.3f dB\n", ...
    name, ...
    metrics.(name).finalNMSD_dB, ...
    metrics.(name).finalMSD_dB, ...
    metrics.(name).finalMSE_dB);

end

%% Plot styles
numEq      = numel(equalizerNames);
colors     = turbo(numEq + 2);
colors     = colors(2:end-1, :);        % skip near-black and near-white ends
lineStyles = ["-", "--", "-.", ":"];
markers    = ["o", "s", "d", "^", "v", "p", "h", ">", "<"];
smoothWin  = 200;                       % moving average (in linear scale) for curves

%% Learning curves (smoothed)
curveNames = ["nmsd", "msd", "mse"];
curveLabel = ["NMSD", "MSD", "MSE"];

for c = 1:numel(curveNames)

    figure;
    hold on;

    for k = 1:numEq

        name = equalizerNames(k);
        lin  = 10.^(metrics.(name).(curveNames(c) + "_dB") / 10);
        lin(~isfinite(lin)) = NaN;

        plot(10 * log10(movmean(lin, smoothWin, 'omitnan')), ...
            'Color', colors(k, :), ...
            'LineStyle', lineStyles(mod(k - 1, numel(lineStyles)) + 1), ...
            'LineWidth', 1.4);

    end

    grid on;
    xlabel('Iteration');
    ylabel(curveLabel(c) + " (dB)");
    title(sprintf('%s - %s (moving average over %d samples)', ...
          data.name, curveLabel(c), smoothWin));

    legend(equalizerNames, ...
        'Interpreter', 'none', ...
        'Location', 'eastoutside');

    hold off;

end

%% Steady-state NMSD (one point per equalizer)

finalNMSD_dB = zeros(numEq, 1);
for k = 1:numEq
    finalNMSD_dB(k) = metrics.(equalizerNames(k)).finalNMSD_dB;
end

figure;
hold on;

for k = 1:numEq
    plot(k, finalNMSD_dB(k), ...
        'Marker', markers(mod(k - 1, numel(markers)) + 1), ...
        'MarkerSize', 9, ...
        'MarkerFaceColor', colors(k, :), ...
        'MarkerEdgeColor', colors(k, :), ...
        'LineStyle', 'none');
    text(k, finalNMSD_dB(k), sprintf('  %.1f', finalNMSD_dB(k)), ...
        'VerticalAlignment', 'middle', 'FontSize', 9);
end

xlim([0.5, numEq + 0.5]);
xticks(1:numEq);
xticklabels(equalizerNames);
set(gca, 'TickLabelInterpreter', 'none');
xtickangle(30);
grid on;
ylabel('Steady-state NMSD (dB)');
title(sprintf('%s - steady-state NMSD (last 10%% of iterations)', data.name));

hold off;

%% True channel evolution

figure;

plot(data.h.');

grid on;
xlabel('Iteration');
ylabel('True channel coefficients');
title('True Channel Evolution');