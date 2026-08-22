function data = generateReceivedData(d, cfg)
%GENERATERECEIVEDDATA Generate received data from composable factors.
%
% The input sequence d is supplied by the caller. This function owns only the
% channel trajectory, receiver noise, and received signal x.

h = generateChannelTrajectory(cfg);

xClean = applyTimeVaryingFIR(d, h);
x = addNoise(xClean, cfg.noise);

data.d = d;
data.x = x;
data.h = h;
data.changeIdx = [];
data.name = describeScenario(cfg);
data.cfg = cfg;
end

function name = describeScenario(cfg)
initializerName = func2str(cfg.channel.initializer);

if isempty(cfg.channel.factors)
    factorText = "static";
else
    factorNames = strings(1, numel(cfg.channel.factors));
    for k = 1:numel(cfg.channel.factors)
        factorNames(k) = string(func2str(cfg.channel.factors{k}));
    end
    factorText = strjoin(factorNames, "+");
end

name = initializerName + "+" + factorText;
end
