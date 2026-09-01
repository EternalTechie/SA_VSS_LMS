%% SA-VSS-LMS Adaptive Filtering Simulation
% Main entry point.
% Change the configured signal, channel, noise, and adaptive filter factors
% here to run a different combination under the same evaluation framework.

clear; clc; close all;

addpath(genpath(pwd));

%% Select adaptive filter

equalizerNames = ["nlms", "vss_lms", "za_lms", "rza_lms"];

%% Global simulation parameters

cfg.N          = 60000;     % number of samples
cfg.M          = 32;        % adaptive filter length
cfg.SNR_dB     = 30;        % received-signal SNR
cfg.delay      = 0;         % desired-signal alignment

% Discrete-time input/training signal
cfg.signal.generator = @whiteGaussianSignal;

% Channel model
cfg.channel.initializer = @sparseChannel;
cfg.channel.M           = cfg.M;
cfg.channel.activeTaps  = [3 9 17];
cfg.channel.gains       = [0.8 -0.5 0.3];
cfg.channel.normalize   = true;
cfg.channel.slowDrift.stdPerSample = 1e-4;
cfg.channel.factors     = {@slowDriftChannel};        % Example: {@slowDriftChannel}

% Receiver noise model
cfg.noise.apply  = @awgnNoise;
cfg.noise.SNR_dB = cfg.SNR_dB;

% Reproducibility
cfg.seed        = 1;

%% Algorithm parameters

% LMS
alg.lms_c.mu    = 0.05;

% NLMS
alg.nlms.mu    = 0.5;
alg.nlms.delta = 1e-8;

% VSS-LMS
alg.vss_lms.mu0   = 0.01;
alg.vss_lms.alpha = 0.97;
alg.vss_lms.gamma = 1e-3;
alg.vss_lms.muMin = 1e-5;
alg.vss_lms.muMax = 0.1;

% ZA-LMS
alg.za_lms.mu  = 0.009;
alg.za_lms.rho = 4e-5;

% RZA-LMS
alg.rza_lms.mu    = 0.009;
alg.rza_lms.rho   = 4e-5;
alg.rza_lms.alpha = 10;

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