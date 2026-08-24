function x = addNoise(xClean, noiseCfg)
%ADDNOISE Apply the configured receiver noise model.

if ~isfield(noiseCfg, "apply")
    error("cfg.noise.apply must be specified.");
end

x = noiseCfg.apply(xClean, noiseCfg);
end
