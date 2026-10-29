%% Monte Carlo Study
% Averages NMSD over many random sparse channels with random sparsity.
% Calls main.m once to initialize cfg and alg, then runs every adaptive
% filter on the same random channels, inputs, and noise.

clear; clc; close all;

%% Make sure rho tuning results exist
% Runs rho_sweep.m first if its results file is missing.
% Set forceRetune = true to rerun it even if the file exists.
forceRetune = false;

if forceRetune || ~isfile("rho_sweep_results.mat")
    fprintf("Running rho_sweep.m first...\n\n");
    run('rho_sweep.m');
    clear; close all;
    fprintf("\nrho_sweep.m finished. Starting Monte Carlo study...\n\n");
end

%% Initialize the main simulation
% If main.m contains the SAVSS_INIT_ONLY check, it returns right after
% defining cfg and alg. Otherwise it runs fully once, which also works.
setappdata(0, 'SAVSS_INIT_ONLY', true);
run('main.m');
if isappdata(0, 'SAVSS_INIT_ONLY')
    rmappdata(0, 'SAVSS_INIT_ONLY');
end
close all;

%% Algorithms to compare
mcNames = ["lms_c", "nlms", "pnlms", "ipnlms", ...
           "za_lms", "rza_lms", "l0_lms", "pza_lms", "cza_lms", ...
           "mvc", "andrew_sine", "hcaf", "lhcaf", "gsrs"];

%% Monte Carlo parameters
mc.N          = 20000;      % iterations per channel
mc.Nch        = 400;        % number of random channels
mc.Kmin       = 2;          % minimum number of active taps
mc.maxDensity = 0.25;       % maximum fraction of active taps
mc.hMin       = 0.2;        % minimum active tap magnitude (threshold)
mc.hMax       = 1.0;        % maximum active tap magnitude
mc.tailFrac   = 0.1;        % final NMSD = mean over last 10% of iterations
mc.seedBase   = 0;          % rho_sweep.m uses a different base (no overlap)
mc.divergeNMSD = 1;         % final NMSD above 0 dB (or non-finite) = diverged

% Convergence-speed metrics (computed on the averaged NMSD curve)
mc.convMargin_dB = 3;       % T_conv: first iteration within 3 dB of own floor
mc.target_dB     = -30;     % T_target: first iteration below -30 dB
mc.smoothWin     = 50;      % moving-average window before thresholding

% Use rho values tuned by rho_sweep.m if available
mc.useTunedRho = true;
mc.rhoFile     = "rho_sweep_results.mat";

%% Apply simulation settings
cfg.N = mc.N;
cfg.channel.normalize = false;    % NMSD is already normalized by ||h||^2
cfg.channel.factors   = {};       % static channel; e.g. {@slowDriftChannel}

mc.Kmax = max(mc.Kmin, floor(mc.maxDensity * cfg.M));
Kvals   = mc.Kmin:mc.Kmax;

if alg.cza_lms.a >= mc.hMin
    warning("CZA cutoff a = %.3g is not below hMin = %.3g.", ...
            alg.cza_lms.a, mc.hMin);
end

if mc.useTunedRho && isfile(mc.rhoFile)
    S = load(mc.rhoFile);
    tunedNames = fieldnames(S.bestRho);

    % Check that rho was tuned under the same settings
    mismatch = strings(0);
    if ~isfield(S, "rs") || ~isfield(S, "alg") || ~isfield(S, "cfg")
        mismatch(end+1) = "results file is from an older rho_sweep.m";
    else
        for f = ["N", "Kmin", "maxDensity", "hMin", "hMax"]
            if ~isfield(S.rs, f) || ~isequal(S.rs.(f), mc.(f))
                mismatch(end+1) = "channel setting " + f;
            end
        end
        for f = ["M", "SNR_dB"]
            if ~isequal(S.cfg.(f), cfg.(f))
                mismatch(end+1) = "cfg." + f;
            end
        end
        for i = 1:numel(tunedNames)
            name = tunedNames{i};
            if ~isfield(S.alg, name) || ~isfield(alg, name) || ...
               ~isequal(rmfield(S.alg.(name), "rho"), rmfield(alg.(name), "rho"))
                mismatch(end+1) = "parameters of " + name;
            end
        end
    end

    for i = 1:numel(mcNames)
        name = char(mcNames(i));
        if isfield(alg, name) && isfield(alg.(name), "rho") && ...
           ~isfield(S.bestRho, name)
            mismatch(end+1) = "no tuned rho for " + name;
        end
    end

    if ~isempty(mismatch)
        warning("Tuned rho may not match current settings. Differs: %s. " + ...
                "Align rho_sweep.m settings and set forceRetune = true.", ...
                strjoin(mismatch, ", "));
    end

    for i = 1:numel(tunedNames)
        if isfield(alg, tunedNames{i})
            alg.(tunedNames{i}).rho = S.bestRho.(tunedNames{i});
            fprintf("Using tuned rho for %s: %.3g\n", ...
                    tunedNames{i}, S.bestRho.(tunedNames{i}));
        end
    end
else
    fprintf("Using rho values from main.m\n");
end

%% Storage
nAlg = numel(mcNames);
nK   = numel(Kvals);

sumNMSD   = zeros(mc.N, nAlg);         % sum of linear NMSD curves, all channels
sumNMSDK  = zeros(mc.N, nAlg, nK);     % same, per K
countK    = zeros(nK, 1);
validCnt  = zeros(1, nAlg);            % non-diverged channels per algorithm
validCntK = zeros(nAlg, nK);           % same, per K
finalNMSD = nan(mc.Nch, nAlg);         % per-channel final NMSD (linear)
diverged  = false(mc.Nch, nAlg);       % true if the run diverged
Kch       = nan(mc.Nch, 1);

tailIdx = floor((1 - mc.tailFrac) * mc.N) + 1 : mc.N;

%% Run Monte Carlo
tStart = tic;

for ch = 1:mc.Nch

    % Same channel, input, and noise for every algorithm
    rng(mc.seedBase + ch);

    [taps, gains, K] = drawSparseChannel(cfg.M, mc);
    cfg.channel.activeTaps = taps;
    cfg.channel.gains      = gains;
    Kch(ch) = K;
    kIdx = K - mc.Kmin + 1;
    countK(kIdx) = countK(kIdx) + 1;

    d    = cfg.signal.generator(cfg.N, cfg.signal);
    data = generateReceivedData(d, cfg);
    hPow = sum(data.h.^2, 1);

    for k = 1:nAlg
        res  = runEqualizer(mcNames(k), data, cfg, alg);
        nmsd = (sum((res.w - data.h).^2, 1) ./ hPow).';

        final = mean(nmsd(tailIdx));

        % Diverged runs are counted but kept out of the averages
        if any(~isfinite(nmsd)) || final > mc.divergeNMSD
            diverged(ch, k) = true;
            continue;
        end

        sumNMSD(:, k)        = sumNMSD(:, k) + nmsd;
        sumNMSDK(:, k, kIdx) = sumNMSDK(:, k, kIdx) + nmsd;
        validCnt(k)          = validCnt(k) + 1;
        validCntK(k, kIdx)   = validCntK(k, kIdx) + 1;
        finalNMSD(ch, k)     = final;
    end

    if mod(ch, 10) == 0
        fprintf("Channel %d/%d | %.0f s elapsed\n", ch, mc.Nch, toc(tStart));
    end

end

%% Averaged curves
avgNMSD  = sumNMSD ./ max(validCnt, 1);
avgNMSDK = sumNMSDK ./ reshape(max(validCntK, 1), 1, nAlg, nK);

%% Final NMSD and convergence metrics
finalMean   = mean(finalNMSD, 1, 'omitnan').';
finalSEM    = std(finalNMSD, 0, 1, 'omitnan').' ./ sqrt(max(validCnt, 1)).';
finalAll_dB = 10 * log10(finalNMSD);
Diverged    = sum(diverged, 1).';

FinalNMSD_dB = 10 * log10(finalMean);
CI95_low_dB  = 10 * log10(max(finalMean - 1.96 * finalSEM, eps));
CI95_high_dB = 10 * log10(finalMean + 1.96 * finalSEM);
Std_dB       = std(finalAll_dB, 0, 1, 'omitnan').';

T_conv   = nan(nAlg, 1);
T_target = nan(nAlg, 1);
for k = 1:nAlg
    [T_conv(k), T_target(k)] = convergenceMetrics(avgNMSD(:, k), mc);
end

summaryTable = table(mcNames(:), FinalNMSD_dB, CI95_low_dB, CI95_high_dB, ...
                     Std_dB, T_target, T_conv, Diverged, ...
    'VariableNames', {'Algorithm', 'FinalNMSD_dB', 'CI95_low_dB', ...
                      'CI95_high_dB', 'Std_dB', 'T_target', 'T_conv', 'Diverged'});
summaryTable = sortrows(summaryTable, 'FinalNMSD_dB');

fprintf("\n=== All channels (%d channels, K = %d to %d, |h| in [%.2g, %.2g]) ===\n", ...
        mc.Nch, mc.Kmin, mc.Kmax, mc.hMin, mc.hMax);
fprintf("FinalNMSD_dB : mean over channels (linear), then dB\n");
fprintf("CI95         : 95%% confidence interval of that mean\n");
fprintf("Std_dB       : spread of per-channel final NMSD\n");
fprintf("T_target     : iterations to reach %g dB\n", mc.target_dB);
fprintf("T_conv       : iterations to reach within %g dB of own floor\n", ...
        mc.convMargin_dB);
fprintf("Diverged     : channels where the run diverged (excluded from averages)\n\n");
disp(summaryTable);

%% Final NMSD per sparsity level
perK_dB = nan(nAlg, nK);
for kk = 1:nK
    sel = Kch == Kvals(kk);
    if any(sel)
        perK_dB(:, kk) = 10 * log10(mean(finalNMSD(sel, :), 1, 'omitnan')).';
    end
end

perKTable = array2table(perK_dB, ...
    'RowNames', cellstr(mcNames), ...
    'VariableNames', cellstr("K_" + string(Kvals)));

fprintf("\n=== Final NMSD (dB) per sparsity level ===\n");
fprintf("Channels per K: %s\n\n", mat2str(countK.'));
disp(perKTable);

%% Plot styles
colors     = turbo(nAlg + 2);
colors     = colors(2:end-1, :);        % skip near-black and near-white ends
lineStyles = ["-", "--", "-.", ":"];
plotIdx    = unique([1:10:mc.N, mc.N]);

%% Averaged NMSD across all channels
figure;
hold on;
for k = 1:nAlg
    plot(plotIdx, 10 * log10(avgNMSD(plotIdx, k)), ...
        'Color', colors(k, :), ...
        'LineStyle', lineStyles(mod(k - 1, numel(lineStyles)) + 1), ...
        'LineWidth', 1.4);
end
yline(mc.target_dB, ':', 'Target', 'HandleVisibility', 'off');
grid on;
xlabel('Iteration');
ylabel('Average NMSD (dB)');
title(sprintf('Average NMSD over %d random sparse channels', mc.Nch));
legend(mcNames, 'Interpreter', 'none', 'Location', 'eastoutside');
hold off;

%% Averaged NMSD per sparsity level
figure;
tiledlayout('flow');
for kk = 1:nK
    if countK(kk) == 0
        continue;
    end
    nexttile;
    hold on;
    for k = 1:nAlg
        plot(plotIdx, 10 * log10(avgNMSDK(plotIdx, k, kk)), ...
            'Color', colors(k, :), ...
            'LineStyle', lineStyles(mod(k - 1, numel(lineStyles)) + 1), ...
            'LineWidth', 1.2);
    end
    grid on;
    xlabel('Iteration');
    ylabel('Average NMSD (dB)');
    title(sprintf('K = %d (%d channels)', Kvals(kk), countK(kk)));
    hold off;
end
lg = legend(mcNames, 'Interpreter', 'none');
lg.Layout.Tile = 'east';

%% Final NMSD vs sparsity level
figure;
hold on;
for k = 1:nAlg
    plot(Kvals, perK_dB(k, :), '-o', ...
        'Color', colors(k, :), ...
        'LineStyle', lineStyles(mod(k - 1, numel(lineStyles)) + 1), ...
        'LineWidth', 1.4);
end
grid on;
xlabel('Number of active taps K');
ylabel('Final NMSD (dB)');
title('Final NMSD vs. sparsity level');
legend(mcNames, 'Interpreter', 'none', 'Location', 'eastoutside');
hold off;

%% Final NMSD with confidence intervals
figure;
hold on;
FinalNMSD_dB_unsorted = 10 * log10(finalMean);
errLow  = FinalNMSD_dB_unsorted - 10 * log10(max(finalMean - 1.96 * finalSEM, eps));
errHigh = 10 * log10(finalMean + 1.96 * finalSEM) - FinalNMSD_dB_unsorted;
errorbar(1:nAlg, FinalNMSD_dB_unsorted, errLow, errHigh, 'o', 'LineWidth', 1.4);
xticks(1:nAlg);
xticklabels(mcNames);
set(gca, 'TickLabelInterpreter', 'none');
grid on;
ylabel('Final NMSD (dB)');
title('Final NMSD with 95% confidence intervals');
hold off;

%% Convergence speed
figure;
barh(categorical(mcNames, mcNames), T_target);
set(gca, 'TickLabelInterpreter', 'none');
grid on;
xlabel(sprintf('Iterations to reach %g dB', mc.target_dB));
title('Convergence speed (averaged curve)');

%% Save results
save("monte_carlo_results.mat", "mc", "mcNames", "avgNMSD", "avgNMSDK", ...
     "finalNMSD", "diverged", "Kch", "summaryTable", "perKTable", "alg", "cfg");
fprintf("\nSaved results to monte_carlo_results.mat\n");

%% Local functions

function [taps, gains, K] = drawSparseChannel(M, mc)
%DRAWSPARSECHANNEL Random sparse channel: K active taps at random positions,
% each with magnitude drawn directly from [mc.hMin, mc.hMax] and random sign.

K     = randi([mc.Kmin, mc.Kmax]);
taps  = sort(randperm(M, K));
mags  = mc.hMin + (mc.hMax - mc.hMin) * rand(K, 1);
signs = 2 * (rand(K, 1) >= 0.5) - 1;
gains = mags .* signs;
end

function [tConv, tTarget] = convergenceMetrics(avgLin, mc)
%CONVERGENCEMETRICS Convergence times from an averaged linear NMSD curve.
%   tConv   : first iteration within mc.convMargin_dB of the curve's own floor
%   tTarget : first iteration below mc.target_dB (NaN if never reached)

N        = numel(avgLin);
tail     = floor((1 - mc.tailFrac) * N) + 1 : N;
floor_dB = 10 * log10(mean(avgLin(tail)));

smooth_dB = movmean(10 * log10(avgLin), mc.smoothWin);

tConv = find(smooth_dB <= floor_dB + mc.convMargin_dB, 1);
if isempty(tConv)
    tConv = NaN;
end

tTarget = find(smooth_dB <= mc.target_dB, 1);
if isempty(tTarget)
    tTarget = NaN;
end
end