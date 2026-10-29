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
    case "l0_lms"
        result = l0_lms(data, cfg, alg.l0_lms);
    case "pza_lms"
        result = pza_lms(data, cfg, alg.pza_lms);
    case "cza_lms"
        result = cza_lms(data, cfg, alg.cza_lms);
    case "scza_lms"
        result = scza_lms(data, cfg, alg.scza_lms);

    % Proportionate algorithms
    case "pnlms"
        result = pnlms(data, cfg, alg.pnlms);
    case "ipnlms"
        result = ipnlms(data, cfg, alg.ipnlms);

    % Robust adaptive filtering algorithms
    case "mvc"
        result = mvc(data, cfg, alg.mvc);
    case "andrew_sine"
        result = andrew_sine(data, cfg, alg.andrew_sine);
    case "hcaf"
        result = hcaf(data, cfg, alg.hcaf);
    case "lhcaf"
        result = lhcaf(data, cfg, alg.lhcaf);
    case "gsrs"
        result = gsrs(data, cfg, alg.gsrs);
    case "sine_cza"
        result = sine_cza(data, cfg, alg.sine_cza);
    case "jlhcaf"
        result = jlhcaf(data, cfg, alg.jlhcaf);
        
    otherwise
        error("Unknown equalizer: %s", name);
end
end