function h = generateChannelTrajectory(cfg)
%GENERATECHANNELTRAJECTORY Build h(:,n) from an initializer and factors.

if ~isfield(cfg, "channel") || ~isfield(cfg.channel, "initializer")
    error("cfg.channel.initializer must be specified.");
end

h0 = cfg.channel.initializer(cfg.channel);
h = repmat(h0, 1, cfg.N);

if ~isfield(cfg.channel, "factors") || isempty(cfg.channel.factors)
    return;
end

for k = 1:numel(cfg.channel.factors)
    factor = cfg.channel.factors{k};
    h = factor(h, cfg);
end
end
