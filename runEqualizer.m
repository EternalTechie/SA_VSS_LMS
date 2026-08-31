function result = runEqualizer(name, data, cfg, alg)
%RUNEQUALIZER Dispatch to the selected adaptive equalizer.

switch lower(string(name))
    case "lms_c"
        result = lms_c(data, cfg, alg.lms_c);
    case "nlms"
        result = nlms(data, cfg, alg.nlms);
    case "vss_lms"
        result = vss_lms(data, cfg, alg.vss_lms);
    case "za_lms"
        result = za_lms(data, cfg, alg.za_lms);
    case "rza_lms"
        result = rza_lms(data, cfg, alg.rza_lms);
    otherwise
        error("Unknown equalizer: %s", name);
end
end
