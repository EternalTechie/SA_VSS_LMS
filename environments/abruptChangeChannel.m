function h = abruptChangeChannel(h, cfg)
%ABRUPTCHANGECHANNEL Overwrite channel trajectory after an abrupt change point.

N = size(h, 2);
M = size(h, 1);

if isfield(cfg.channel, "abrupt") && isfield(cfg.channel.abrupt, "changeIdx")
    changeIdx = cfg.channel.abrupt.changeIdx;
else
    changeIdx = floor(N / 2);
end

if isfield(cfg.channel, "abrupt") && isfield(cfg.channel.abrupt, "activeTaps")
    activeTaps = cfg.channel.abrupt.activeTaps;
    gains = cfg.channel.abrupt.gains;
else
    activeTaps = [5 13 27];
    gains = [0.6 0.5 -0.7];
end

h2 = zeros(M, 1);
h2(activeTaps) = gains;
if isfield(cfg.channel, "normalize") && cfg.channel.normalize
    h2 = h2 / norm(h2);
end

h(:, changeIdx+1:end) = repmat(h2, 1, N - changeIdx);
end
