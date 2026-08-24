function h0 = sparseChannel(channelCfg)
%SPARSECHANNEL Create an initial sparse finite impulse response channel.

M = channelCfg.M;
h0 = zeros(M, 1);

h0(channelCfg.activeTaps) = channelCfg.gains;

if isfield(channelCfg, "normalize") && channelCfg.normalize
    hNorm = norm(h0);
    if hNorm > 0
        h0 = h0 / hNorm;
    end
end
end
