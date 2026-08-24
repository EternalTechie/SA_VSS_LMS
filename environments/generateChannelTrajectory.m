function [h, info] = generateChannelTrajectory(cfg)
%GENERATECHANNELTRAJECTORY Build h(:,n) from an initializer and factors.

if ~isfield(cfg, "channel") || ~isfield(cfg.channel, "initializer")
    error("cfg.channel.initializer must be specified.");
end

h0 = cfg.channel.initializer(cfg.channel);
h = repmat(h0, 1, cfg.N);

info.initializer = func2str(cfg.channel.initializer);
info.factors = strings(1, 0);
info.factorDeltaNorm = [];
info.changeIdx = [];

if ~isfield(cfg.channel, "factors") || isempty(cfg.channel.factors)
    return;
end

for k = 1:numel(cfg.channel.factors)
    factor = cfg.channel.factors{k};
    factorName = string(func2str(factor));

    if factorName == "applyTimeVaryingFIR"
        error(["applyTimeVaryingFIR is a channel renderer, not a channel factor. " + ...
            "Use cfg.channel.renderer = @applyTimeVaryingFIR instead."]);
    end

    hBefore = h;
    h = factor(h, cfg);

    info.factors(end+1) = factorName;
    info.factorDeltaNorm(end+1) = norm(h(:) - hBefore(:));
end

if isfield(cfg.channel, "abrupt") && isfield(cfg.channel.abrupt, "changeIdx")
    info.changeIdx = cfg.channel.abrupt.changeIdx;
end
end
