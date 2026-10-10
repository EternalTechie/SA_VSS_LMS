%% CZA-LMS Ablation Study on the Cutoff Parameter a
% Sweeps the cutoff a for one or more channels whose taps are set HERE,
% to show how the influence of a depends on the true tap magnitudes.
% Calls main.m once to initialize cfg and alg; the channel in main.m is
% not used.

clear; clc; close all;

%% Initialize the main simulation
% If main.m contains the SAVSS_INIT_ONLY check, it returns right after
% defining cfg and alg. Otherwise it runs fully once, which also works.
setappdata(0, 'SAVSS_INIT_ONLY', true);
run('main.m');
if isappdata(0, 'SAVSS_INIT_ONLY')
    rmappdata(0, 'SAVSS_INIT_ONLY');
end
close all;

%% Channels to test
% Set tap positions (1..M) and gains for each channel.
% Add or remove entries to compare different tap magnitudes.
% label is optional; leave it as "" to label by the smallest tap.
chans = struct('taps', {}, 'gains', {}, 'label', {});

chans(end+1) = struct('taps', [2 4 9], 'gains', [0.2 0.8 0.7], 'label', "");
chans(end+1) = struct('taps', [2 4 9], 'gains', [0.4 0.8 0.7], 'label', "");
chans(end+1) = struct('taps', [2 4 9], 'gains', [0.8 0.8 0.7], 'label', "");
chans(end+1) = struct('taps', [2 4 9], 'gains', [1.2 0.8 0.7], 'label', "");

%% Ablation parameters
abl.aValues   = [0.01 0.02 0.03 0.05 0.07 0.1 0.15 0.2 0.3 0.5 0.7 1 1.5 2 3];
abl.N         = 20000;      % iterations per run
abl.R         = 20;         % independent input/noise realizations per channel
abl.normalize = false;      % true: scale each channel to unit norm
abl.tailFrac  = 0.1;        % final NMSD = mean over last 10% of iterations
abl.seedBase  = 7e5;
abl.edge_dB   = 1;          % window edge: a where NMSD is 1 dB worse than best

% Rho choice: "main" uses alg.cza_lms.rho from main.m,
% "tuned" uses the best rho from rho_sweep_results.mat
abl.rhoMode = "main";
abl.rhoFile = "rho_sweep_results.mat";

% Plot averaged NMSD convergence curves for every a (one figure per channel)
abl.plotCurves = false;

%% Apply settings
cfg.N = abl.N;
cfg.channel.normalize = abl.normalize;
cfg.channel.factors   = {};       % static channel

if abl.rhoMode == "tuned"
    S = load(abl.rhoFile, "bestRho");
    alg.cza_lms.rho = S.bestRho.cza_lms;
end
fprintf("CZA-LMS: mu = %.3g, rho = %.3g\n", alg.cza_lms.mu, alg.cza_lms.rho);

nC = numel(chans);
nA = numel(abl.aValues);

minTap = zeros(nC, 1);
for c = 1:nC
    taps  = chans(c).taps;
    gains = chans(c).gains;
    if numel(taps) ~= numel(gains)
        error("Channel %d: taps and gains must have the same length.", c);
    end
    if any(taps < 1 | taps > cfg.M) || numel(unique(taps)) ~= numel(taps)
        error("Channel %d: taps must be unique positions between 1 and %d.", c, cfg.M);
    end

    g = gains(:);
    if abl.normalize
        g = g / norm(g);
    end
    minTap(c) = min(abs(g));

    if chans(c).label == ""
        chans(c).label = sprintf("min |h| = %.3g", minTap(c));
    end
end

if max(abl.aValues) <= max(minTap)
    warning("Largest a (%.3g) does not exceed the largest smallest-tap (%.3g); " + ...
            "the upper edge of the window may not be visible.", ...
            max(abl.aValues), max(minTap));
end

%% Storage
finalNMSD = nan(abl.R, nA, nC);           % per-run final NMSD (linear)
sumCurves = zeros(abl.N, nA, nC);         % sum of NMSD curves (linear)

tailIdx = floor((1 - abl.tailFrac) * abl.N) + 1 : abl.N;

%% Run ablation
tStart = tic;

for c = 1:nC

    cfg.channel.activeTaps = chans(c).taps;
    cfg.channel.gains      = chans(c).gains(:);

    for r = 1:abl.R

        % Same input and noise for every channel and every a
        rng(abl.seedBase + r);

        d    = cfg.signal.generator(cfg.N, cfg.signal);
        data = generateReceivedData(d, cfg);
        hPow = sum(data.h.^2, 1);

        for j = 1:nA
            algCZA   = alg.cza_lms;
            algCZA.a = abl.aValues(j);

            res  = cza_lms(data, cfg, algCZA);
            nmsd = (sum((res.w - data.h).^2, 1) ./ hPow).';

            finalNMSD(r, j, c) = mean(nmsd(tailIdx));
            sumCurves(:, j, c) = sumCurves(:, j, c) + nmsd;
        end

    end

    fprintf("Channel %d/%d (%s) done | %.0f s elapsed\n", ...
            c, nC, chans(c).label, toc(tStart));

end

%% Results
meanFinal = reshape(mean(finalNMSD, 1), nA, nC);
semFinal  = reshape(std(finalNMSD, 0, 1), nA, nC) / sqrt(abl.R);

final_dB = 10 * log10(meanFinal);
errLo_dB = final_dB - 10 * log10(max(meanFinal - 1.96 * semFinal, eps));
errHi_dB = 10 * log10(meanFinal + 1.96 * semFinal) - final_dB;

best_dB = nan(nC, 1);
bestA   = nan(nC, 1);
aLow    = nan(nC, 1);
aHigh   = nan(nC, 1);

for c = 1:nC
    [best_dB(c), iBest] = min(final_dB(:, c));
    bestA(c) = abl.aValues(iBest);
    [aLow(c), aHigh(c)] = windowEdges(abl.aValues, final_dB(:, c), iBest, abl.edge_dB);
end

windowTable = table(string({chans.label}).', minTap, best_dB, bestA, aLow, aHigh, ...
                    aHigh ./ minTap, ...
    'VariableNames', {'Channel', 'MinTap', 'BestNMSD_dB', 'BestA', ...
                      'WindowLow', 'WindowHigh', 'High_over_MinTap'});

fprintf("\n=== Safe window for a (within %g dB of best) ===\n", abl.edge_dB);
fprintf("WindowLow / WindowHigh: a values where NMSD is %g dB worse than best\n", abl.edge_dB);
fprintf("High_over_MinTap: upper edge relative to the smallest true tap\n\n");
disp(windowTable);

fprintf("\n=== Final NMSD (dB) vs a ===\n\n");
disp(array2table(final_dB, ...
    'RowNames', cellstr("a_" + string(abl.aValues)), ...
    'VariableNames', cellstr("ch" + string(1:nC))));

%% Plot styles
colors     = turbo(nC + 2);
colors     = colors(2:end-1, :);        % skip near-black and near-white ends
lineStyles = ["-", "--", "-.", ":"];

%% Final NMSD vs a
figure;
hold on;
for c = 1:nC
    errorbar(abl.aValues, final_dB(:, c), errLo_dB(:, c), errHi_dB(:, c), '-o', ...
        'Color', colors(c, :), ...
        'LineStyle', lineStyles(mod(c - 1, numel(lineStyles)) + 1), ...
        'LineWidth', 1.5, ...
        'DisplayName', chans(c).label);
end
for c = 1:nC
    xline(minTap(c), ':', 'Color', colors(c, :), 'LineWidth', 1.2, ...
          'HandleVisibility', 'off');
end
set(gca, 'XScale', 'log');
grid on;
xlabel('Cutoff parameter a');
ylabel('Final NMSD (dB)');
title('CZA-LMS: final NMSD vs. a (dotted lines = smallest true tap)');
legend('Location', 'best');
hold off;

%% Final NMSD vs a / smallest tap
figure;
hold on;
for c = 1:nC
    plot(abl.aValues / minTap(c), final_dB(:, c) - best_dB(c), '-o', ...
        'Color', colors(c, :), ...
        'LineStyle', lineStyles(mod(c - 1, numel(lineStyles)) + 1), ...
        'LineWidth', 1.5, ...
        'DisplayName', chans(c).label);
end
xline(1, 'k:', 'a = smallest tap', 'HandleVisibility', 'off');
yline(abl.edge_dB, 'k--', sprintf('+%g dB', abl.edge_dB), 'HandleVisibility', 'off');
set(gca, 'XScale', 'log');
grid on;
xlabel('a / smallest true tap');
ylabel('NMSD relative to best (dB)');
title('If the upper knees line up, the window scales with the smallest tap');
legend('Location', 'best');
hold off;

%% Upper window edge vs smallest tap
if nC > 1
    figure;
    hold on;
    plot(minTap, aHigh, 'o', 'MarkerSize', 8, 'LineWidth', 1.5, ...
         'DisplayName', 'measured upper edge');
    lims = [0, max([minTap; aHigh], [], 'omitnan') * 1.1];
    plot(lims, lims, 'k:', 'DisplayName', 'a = smallest tap');
    valid = isfinite(aHigh);
    if nnz(valid) >= 2
        p = polyfit(minTap(valid), aHigh(valid), 1);
        plot(lims, polyval(p, lims), '--', 'LineWidth', 1.2, ...
             'DisplayName', sprintf('fit: %.2f x %+.3f', p(1), p(2)));
    end
    grid on;
    xlabel('Smallest true tap |h|_{min}');
    ylabel(sprintf('Upper edge of a window (+%g dB)', abl.edge_dB));
    title('Where a starts to hurt vs. smallest true tap');
    legend('Location', 'northwest');
    hold off;
end

%% Averaged convergence curves for every a
if abl.plotCurves
    plotIdx = unique([1:10:abl.N, abl.N]);
    aColors = turbo(nA + 2);
    aColors = aColors(2:end-1, :);
    for c = 1:nC
        figure;
        hold on;
        for j = 1:nA
            plot(plotIdx, 10 * log10(sumCurves(plotIdx, j, c) / abl.R), ...
                'Color', aColors(j, :), 'LineWidth', 1.2, ...
                'DisplayName', sprintf('a = %.3g', abl.aValues(j)));
        end
        grid on;
        xlabel('Iteration');
        ylabel('Average NMSD (dB)');
        title(sprintf('CZA-LMS convergence for each a (%s)', chans(c).label));
        legend('Location', 'eastoutside');
        hold off;
    end
end

%% Save results
save("ablation_results.mat", "abl", "chans", "minTap", "finalNMSD", ...
     "final_dB", "windowTable", "alg", "cfg");
fprintf("\nSaved results to ablation_results.mat\n");

%% Local functions

function [aLow, aHigh] = windowEdges(aValues, final_dB, iBest, edge_dB)
%WINDOWEDGES a values on each side of the best a where final NMSD first
% becomes edge_dB worse than the best. Interpolated on a log-a scale.
% NaN if the curve never gets that much worse within the tested range.

level = final_dB(iBest) + edge_dB;
la    = log10(aValues(:));

aHigh = NaN;
for j = iBest + 1 : numel(aValues)
    if final_dB(j) >= level
        t = (level - final_dB(j - 1)) / (final_dB(j) - final_dB(j - 1));
        aHigh = 10^(la(j - 1) + t * (la(j) - la(j - 1)));
        break;
    end
end

aLow = NaN;
for j = iBest - 1 : -1 : 1
    if final_dB(j) >= level
        t = (level - final_dB(j + 1)) / (final_dB(j) - final_dB(j + 1));
        aLow = 10^(la(j + 1) + t * (la(j) - la(j + 1)));
        break;
    end
end
end