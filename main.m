%% SA-VSS-LMS Adaptive Filtering Simulation
% Main entry point.
% Change the configured signal, channel, noise, and adaptive filter factors
% here to run a different combination under the same evaluation framework.

clear; clc; close all;

addpath(genpath(pwd));

%% Select adaptive filter

equalizerName = "nlms";

%% Global simulation parameters

cfg.N          = 20000;     % number of samples
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

alg.mu         = 0.5;
alg.delta      = 1e-8;

%% Generate input signal and received data

rng(cfg.seed);

d = generateInputSignal(cfg.N, cfg);
data = generateReceivedData(d, cfg);

%% Run adaptive equalizer

result = runEqualizer(equalizerName, data, cfg, alg);

%% Evaluate

metrics = evaluateResult(data, result, cfg);

%% Display basic results

fprintf("\nScenario    : %s\n", data.name);
fprintf("Equalizer   : %s\n", equalizerName);
fprintf("Final NMSD  : %.3f dB\n", metrics.finalNMSD_dB);
fprintf("Final error : %.3e\n", metrics.finalErrorPower);

%% Basic plots

figure;
plot(metrics.nmsd_dB, 'LineWidth', 1.2);
grid on;
xlabel('Iteration');
ylabel('NMSD (dB)');
title(sprintf('%s - %s', data.name, equalizerName));

figure;
plot(data.h.');
grid on;
xlabel('Iteration');
ylabel('True channel coefficients');
title('True channel evolution');

figure;
plot(result.w.');
grid on;
xlabel('Iteration');
ylabel('Estimated filter coefficients');
title('Estimated coefficients');
