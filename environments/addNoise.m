function x = addNoise(xClean, noiseCfg)
%ADDNOISE Apply the configured receiver noise model, then any noise factors.
%
% cfg.noise.apply   : base noise model, e.g. @awgnNoise
% cfg.noise.factors : optional cell array of noise factors applied in order
%                     on top of the base noise, e.g. {@impulseNoise}.
%                     Each factor has the form x = factor(x, xClean, noiseCfg).

if ~isfield(noiseCfg, "apply")
    error("cfg.noise.apply must be specified.");
end

x = noiseCfg.apply(xClean, noiseCfg);

if isfield(noiseCfg, "factors") && ~isempty(noiseCfg.factors)
    for k = 1:numel(noiseCfg.factors)
        x = noiseCfg.factors{k}(x, xClean, noiseCfg);
    end
end
end
