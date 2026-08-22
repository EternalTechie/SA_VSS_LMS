function metrics = evaluateResult(data, result, cfg)
%EVALUATERESULT Calculate common evaluation quantities.
%
% This function deliberately does not decide what "good" means. It returns
% quantities that can later be used to define convergence, sparse recovery,
% tracking, and complexity metrics.

N = cfg.N;

nmsd = nan(N,1);

for n = 1:N
    htrue = data.h(:,n);
    hest  = result.w(:,n);

    denom = htrue' * htrue;

    if denom > 0
        nmsd(n) = ((hest-htrue)'*(hest-htrue)) / denom;
    end
end

% Ignore the initial undefined/zero-update portion when displaying results.
valid = isfinite(nmsd) & nmsd > 0;

nmsd_dB = nan(size(nmsd));
nmsd_dB(valid) = 10*log10(nmsd(valid));

% Final value: average over the last 10% of valid iterations.
idx = find(valid);
tailStart = max(1, floor(0.9*numel(idx)));
tailIdx = idx(tailStart:end);

metrics.nmsd = nmsd;
metrics.nmsd_dB = nmsd_dB;
metrics.finalNMSD = mean(nmsd(tailIdx), 'omitnan');
metrics.finalNMSD_dB = 10*log10(metrics.finalNMSD);

metrics.errorPower = result.e.^2;
metrics.finalErrorPower = mean(metrics.errorPower(tailIdx), 'omitnan');
end
