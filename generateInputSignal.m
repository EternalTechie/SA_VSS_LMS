function d = generateInputSignal(N, cfg)
%GENERATEINPUTSIGNAL Generate the simulator input/training sequence.
%
% The signal model is selected outside the channel environment. This keeps
% the transmitted/input sequence independent from receiver/channel effects.

if ~isfield(cfg, "signal") || ~isfield(cfg.signal, "generator")
    error("cfg.signal.generator must be specified.");
end

d = cfg.signal.generator(N, cfg.signal);
end
