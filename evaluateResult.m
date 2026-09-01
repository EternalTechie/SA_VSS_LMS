function metrics = evaluateResult(data, result, cfg)
%EVALUATERESULT Calculate common evaluation quantities.

N = cfg.N;

nmsd = nan(N,1);
msd  = nan(N,1);
mse  = nan(N,1);

for n = 1:N

    htrue = data.h(:,n);
    hest  = result.w(:,n);

    % Mean-square deviation
    diff = hest - htrue;
    msd(n) = diff' * diff;

    % Normalized mean-square deviation
    denom = htrue' * htrue;

    if denom > 0
        nmsd(n) = msd(n) / denom;
    end

    % Instantaneous squared error
    mse(n) = result.e(n)^2;

end

% Ignore invalid initial portion
valid = isfinite(nmsd) & nmsd > 0;

nmsd_dB = nan(size(nmsd));
nmsd_dB(valid) = 10*log10(nmsd(valid));

% Final value: average over last 10% of valid iterations
idx = find(valid);

tailStart = max(1, floor(0.9*numel(idx)));
tailIdx = idx(tailStart:end);

% Store trajectories
metrics.msd = msd;
metrics.mse = mse;

metrics.msd_dB = 10*log10(msd);
metrics.mse_dB = 10*log10(mse);

metrics.nmsd = nmsd;
metrics.nmsd_dB = nmsd_dB;

% Final metrics
metrics.finalMSD = mean(msd(tailIdx), 'omitnan');
metrics.finalMSD_dB = 10*log10(metrics.finalMSD);

metrics.finalMSE = mean(mse(tailIdx), 'omitnan');
metrics.finalMSE_dB = 10*log10(metrics.finalMSE);

metrics.finalNMSD = mean(nmsd(tailIdx), 'omitnan');
metrics.finalNMSD_dB = 10*log10(metrics.finalNMSD);

end