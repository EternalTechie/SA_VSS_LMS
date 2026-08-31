function result = runEqualizer(name, data, cfg, alg)
%RUNEQUALIZER Dispatch to the selected adaptive equalizer.

switch lower(string(name))
    case "nlms"
        result = nlms(data, cfg, alg.nlms);
    case "vss_lms"
        result = vss_lms(data, cfg, alg);
    otherwise
        error("Unknown equalizer: %s", name);
end
end
