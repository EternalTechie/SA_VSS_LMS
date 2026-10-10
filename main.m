%% SA-VSS-LMS Adaptive Filtering Simulation
% Main entry point.
% Change the configured signal, channel, noise, and adaptive filter factors
% here to run a different combination under the same evaluation framework.

clear; clc; close all;

addpath(genpath(pwd));

%% Select adaptive filter

equalizerNames = ["nlms", "vss_lms", "za_lms", "rza_lms", ...
                  "pza_lms", "cza_lms", ...
                  "mvc", "andrew_sine", "hcaf", "lhcaf"];

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
cfg.noise.apply  = @awgnNoise;
cfg.noise.SNR_dB = cfg.SNR_dB;

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
% Hyperbolic Cosine Adaptive Filter
alg.hcaf.lambda = 1;
alg.hcaf.mu     = mu_common / alg.hcaf.lambda;       % 0.01, effective step mu*lambda
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

%% NMSD comparison

figure;
hold on;

for k = 1:numel(equalizerNames)

    name = equalizerNames(k);

    plot(metrics.(name).nmsd_dB, 'LineWidth', 1.2);

end

grid on;
xlabel('Iteration');
ylabel('NMSD (dB)');
title(sprintf('%s - NMSD Comparison', data.name));

legend(equalizerNames, ...
    'Interpreter', 'none', ...
    'Location', 'best');

hold off;

%% MSD comparison

figure;
hold on;

for k = 1:numel(equalizerNames)

    name = equalizerNames(k);

    plot(metrics.(name).msd_dB, 'LineWidth', 1.2);

end

grid on;
xlabel('Iteration');
ylabel('MSD (dB)');
title(sprintf('%s - MSD Comparison', data.name));

legend(equalizerNames, ...
    'Interpreter', 'none', ...
    'Location', 'best');

hold off;

%% MSE comparison

figure;
hold on;

for k = 1:numel(equalizerNames)

    name = equalizerNames(k);

    plot(metrics.(name).mse_dB, 'LineWidth', 1.2);

end

grid on;
xlabel('Iteration');
ylabel('MSE (dB)');
title(sprintf('%s - MSE Comparison', data.name));

legend(equalizerNames, ...
    'Interpreter', 'none', ...
    'Location', 'best');

hold off;

%% True channel evolution

figure;

plot(data.h.');

grid on;
xlabel('Iteration');
ylabel('True channel coefficients');
title('True Channel Evolution');

%% Equalizer coefficient evolution

for k = 1:numel(equalizerNames)

    name = equalizerNames(k);

    figure;

    plot(results.(name).w.');

    grid on;
    xlabel('Iteration');
    ylabel('Estimated channel coefficients');
    title(sprintf('%s - Coefficient Evolution', name));

end