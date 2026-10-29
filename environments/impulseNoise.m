function x = impulseNoise(x, xClean, noiseCfg)
%IMPULSENOISE Add impulsive noise on top of the base receiver noise.
%
% Noise factor: applied after cfg.noise.apply (e.g. awgnNoise), the same
% way channel factors like slowDriftChannel are applied after the channel
% initializer. Enable with cfg.noise.factors = {@impulseNoise}.
%
% Bernoulli-Gaussian impulses: each sample is hit with probability prob
% by an extra Gaussian impulse whose variance is powerRatio times the
% clean signal power.
%
% Optional fixed outliers: samples in outlierIdx get an extra spike of
% magnitude outlierAmp with a random sign.
%
% Settings, all optional, in noiseCfg.impulse:
%   prob       : impulse probability per sample          (default 0.01)
%   powerRatio : impulse variance / clean signal power   (default 1000)
%   outlierIdx : sample indices of fixed outliers        (default [])
%   outlierAmp : magnitude of each fixed outlier         (default 100)
%
% Model follows the Bernoulli-Gaussian setup used by Kumar, Bhattacharjee
% and George, IEEE TCAS-II, 2021 (JLHCAF).

if isfield(noiseCfg, "impulse")
    imp = noiseCfg.impulse;
else
    imp = struct();
end

prob       = getField(imp, "prob", 0.01);
powerRatio = getField(imp, "powerRatio", 1000);
outlierIdx = getField(imp, "outlierIdx", []);
outlierAmp = getField(imp, "outlierAmp", 100);

% Bernoulli-Gaussian impulses
signalPower = mean(xClean.^2);
hit = rand(size(x)) < prob;
x(hit) = x(hit) + sqrt(powerRatio * signalPower) * randn(size(x(hit)));

% Fixed outliers
outlierIdx = outlierIdx(outlierIdx >= 1 & outlierIdx <= numel(x));
if ~isempty(outlierIdx)
    signs = 2 * (rand(size(x(outlierIdx))) >= 0.5) - 1;
    x(outlierIdx) = x(outlierIdx) + outlierAmp * signs;
end
end

function value = getField(s, name, default)
%GETFIELD Return s.(name) if present, otherwise the default.
if isfield(s, name)
    value = s.(name);
else
    value = default;
end
end