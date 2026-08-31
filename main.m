%% SA-VSS-LMS Adaptive Filtering Simulation
% Main entry point.
% Change the configured signal, channel, noise, and adaptive filter factors
% here to run a different combination under the same evaluation framework.

clear; clc; close all;

addpath(genpath(pwd));

%% Select adaptive filter

equalizerNames = ["nlms", "vss_lms"];

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
cfg.channel.factors     = {};        % Example: {@slowDriftChannel}

% Receiver noise model
cfg.noise.apply  = @awgnNoise;
cfg.noise.SNR_dB = cfg.SNR_dB;

% Reproducibility
cfg.seed        = 1;

%% Algorithm parameters

alg.nlms.mu         = 0.5;
alg.nlms.delta      = 1e-8;

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

%% Evaluate

metrics = struct();

for k = 1:numel(equalizerNames)

    name = equalizerNames(k);

    metrics.(name) = evaluateResult(data, results.(name), cfg);

end

%% Display basic results

fprintf("\nScenario: %s\n\n", data.name);

for k = 1:numel(equalizerNames)

    name = equalizerNames(k);
    metrics = evaluateResult(data, results.(name), cfg);

    fprintf("%-10s | Final NMSD: %8.3f dB | Final error: %.3e\n", ...
        name, ...
        metrics.finalNMSD_dB, ...
        metrics.finalErrorPower);

end

%% NMSD comparison

figure;
hold on;

for k = 1:numel(equalizerNames)

    name = equalizerNames(k);

    metrics = evaluateResult(data, results.(name), cfg);

    plot(metrics.nmsd_dB, 'LineWidth', 1.2);

end

grid on;
xlabel('Iteration');
ylabel('NMSD (dB)');
title(sprintf('%s - NMSD Comparison', data.name));

legend(equalizerNames, 'Interpreter', 'none', 'Location', 'best');

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
    result = results.(name);

    figure;
    plot(result.w.');

    grid on;
    xlabel('Iteration');
    ylabel('Equalizer coefficients');
    title(sprintf('%s - Coefficient Evolution', name));

end