function h = slowDriftChannel(h, cfg)
%SLOWDRIFTCHANNEL Add random-walk gain drift to initially active taps.

if isfield(cfg.channel, "slowDrift") && isfield(cfg.channel.slowDrift, "stdPerSample")
    stdPerSample = cfg.channel.slowDrift.stdPerSample;
else
    stdPerSample = 1e-4;
end

active = find(abs(h(:,1)) > 0);
N = size(h, 2);

for idx = 1:numel(active)
    tap = active(idx);
    drift = cumsum(stdPerSample * randn(1, N));
    h(tap, :) = h(tap, :) + drift;
end

if isfield(cfg.channel, "normalize") && cfg.channel.normalize
    for n = 1:N
        hNorm = norm(h(:,n));
        if hNorm > 0
            h(:,n) = h(:,n) / hNorm;
        end
    end
end
end
