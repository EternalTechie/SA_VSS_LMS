%% Small Tap Study
% Tests how each zero-attraction shape treats a real tap that is small
% enough to sit inside (or near) the attraction zone.
%
% Each channel has a few large taps plus ONE small tap of magnitude s.
% s is swept from well inside the zone to just outside it. For every s,
% the script measures:
%   - final NMSD of the whole channel
%   - relative bias of the small tap: (E[w_small] - s) / s
%   - how often the small tap is suppressed (estimate below half its value)
% and compares them with the bias predicted from each attraction function.
%
% Thresholding view (support recovery):
%   An ideal threshold filter REMOVES taps far below the cutoff (treated as
%   noise), KEEPS taps around the cutoff (treated as real), and leaves
%   true-zero taps clean. The script scores each method on all three:
%   - Removal score : how much of a far-below-cutoff tap is removed
%   - Keep score    : how much of a near-cutoff tap is kept
%   - Zero-tap energy: leftover energy on true-zero taps (false taps)
%
% Calls main.m once to initialize cfg and alg.

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

%% Algorithms to compare
stNames = ["lms_c", "l0_lms", "pza_lms", "cza_lms"];

%% Study parameters
st.smallVals = [0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10 0.12 0.15 0.20];
st.Nch       = 100;         % channels per small-tap value
st.N         = 10000;       % iterations per channel
st.KlargeMin = 2;           % minimum number of large taps
st.KlargeMax = 5;           % maximum number of large taps
st.largeMin  = 0.4;         % minimum large tap magnitude
st.largeMax  = 1.0;         % maximum large tap magnitude
st.tailFrac  = 0.1;         % steady state = last 10% of iterations
st.seedBase  = 9e5;         % different channels from the other studies
st.suppressFrac = 0.5;      % small tap counted as suppressed below 50% of s

% Thresholding zones, as fractions of the cutoff
st.removeBelow = 0.5;       % s <= 0.5*cutoff : far below, should be removed
st.keepFrom    = 0.7;       % 0.7*cutoff <= s <= cutoff : near cutoff, should be kept

% Rho choice
%   "common": same rho for all sparse algorithms. All attraction functions
%             peak at 1, so this compares the SHAPE only.
%   "tuned" : best rho per algorithm from rho_sweep_results.mat.
st.rhoMode = "tuned";
st.rhoFile = "rho_sweep_results.mat";

%% Apply simulation settings
cfg.N = st.N;
cfg.channel.normalize = false;    % NMSD is already normalized by ||h||^2
cfg.channel.factors   = {};       % static channel

nAlg   = numel(stNames);
nS     = numel(st.smallVals);
sparse = stNames(stNames ~= "lms_c");

switch st.rhoMode
    case "common"
        for name = sparse
            alg.(name).rho = rho_common;
        end
        fprintf("Using common rho = %.3g for all sparse algorithms\n", rho_common);
    case "tuned"
        S = load(st.rhoFile, "bestRho");
        for name = sparse
            alg.(name).rho = S.bestRho.(name);
            fprintf("Using tuned rho for %s: %.3g\n", name, alg.(name).rho);
        end
    otherwise
        error("st.rhoMode must be ""common"" or ""tuned"".");
end

% All three cutoffs should match so only the shape differs
cutoffs = [alg.cza_lms.a, 1 / alg.l0_lms.beta, 1 / (alg.pza_lms.alpha - 1)];
fprintf("Cutoffs: CZA a = %.3g | l0 1/beta = %.3g | PZA 1/(alpha-1) = %.3g\n", cutoffs);
if max(cutoffs) - min(cutoffs) > 1e-12
    warning("Cutoffs differ, so shape and width are both changing.");
end
cutoff = alg.cza_lms.a;

%% Storage
finalNMSD  = nan(st.Nch, nS, nAlg);   % final NMSD (linear)
relBias    = nan(st.Nch, nS, nAlg);   % (E[w_small] - s) / s
suppressed = false(st.Nch, nS, nAlg); % small tap estimate below suppressFrac*s
zeroEnergy = nan(st.Nch, nS, nAlg);   % energy on true-zero taps / ||h||^2

tailIdx = floor((1 - st.tailFrac) * st.N) + 1 : st.N;

%% Run study
tStart = tic;

for sIdx = 1:nS

    s = st.smallVals(sIdx);

    for ch = 1:st.Nch

        % Same channel, input, and noise for every algorithm
        rng(st.seedBase + 1000 * sIdx + ch);

        Klarge = randi([st.KlargeMin, st.KlargeMax]);
        taps   = randperm(cfg.M, Klarge + 1);
        mags   = [s; st.largeMin + (st.largeMax - st.largeMin) * rand(Klarge, 1)];
        signs  = 2 * (rand(Klarge + 1, 1) >= 0.5) - 1;
        gains  = mags .* signs;

        smallPos  = taps(1);
        smallSign = signs(1);
        zeroIdx   = setdiff(1:cfg.M, taps);

        cfg.channel.activeTaps = taps;
        cfg.channel.gains      = gains;

        d    = cfg.signal.generator(cfg.N, cfg.signal);
        data = generateReceivedData(d, cfg);
        hPow = sum(data.h.^2, 1);

        for k = 1:nAlg
            res  = runEqualizer(stNames(k), data, cfg, alg);
            nmsd = sum((res.w - data.h).^2, 1) ./ hPow;

            finalNMSD(ch, sIdx, k) = mean(nmsd(tailIdx));

            wSmall = smallSign * mean(res.w(smallPos, tailIdx));
            relBias(ch, sIdx, k)    = (wSmall - s) / s;
            suppressed(ch, sIdx, k) = wSmall < st.suppressFrac * s;

            zeroEnergy(ch, sIdx, k) = mean(sum(res.w(zeroIdx, tailIdx).^2, 1)) / ...
                                      sum(gains.^2);
        end

    end

    fprintf("s = %.3f done | %.0f s elapsed\n", s, toc(tStart));

end

%% Results
meanFinal = squeeze(mean(finalNMSD, 1));              % nS x nAlg (linear)
semFinal  = squeeze(std(finalNMSD, 0, 1)) / sqrt(st.Nch);

final_dB   = 10 * log10(meanFinal);
finalLo_dB = final_dB - 10 * log10(max(meanFinal - 1.96 * semFinal, eps));
finalHi_dB = 10 * log10(meanFinal + 1.96 * semFinal) - final_dB;

meanBias = squeeze(mean(relBias, 1)) * 100;           % percent
ciBias   = 1.96 * squeeze(std(relBias, 0, 1)) / sqrt(st.Nch) * 100;

suppressPct = squeeze(mean(suppressed, 1)) * 100;

% Predicted steady-state bias: mu*sx2*(s - w) = rho*c(w), evaluated at w ~ s
sx2 = 1;
predBias = zeros(nS, nAlg);
for k = 1:nAlg
    if stNames(k) == "lms_c"
        continue;
    end
    c = attraction(stNames(k), st.smallVals, alg);
    predBias(:, k) = -100 * alg.(stNames(k)).rho * c(:) ./ ...
                     (alg.(stNames(k)).mu * sx2 * st.smallVals(:));
end
predBias = max(predBias, -100);

% Thresholding scores
keepRatio = 1 + meanBias / 100;                       % E[w_small] / s, nS x nAlg

removeZone = st.smallVals <= st.removeBelow * cutoff;
keepZone   = st.smallVals >= st.keepFrom * cutoff & st.smallVals <= cutoff;

RemovalScore = 100 * mean(min(max(1 - keepRatio(removeZone, :), 0), 1), 1).';
KeepScore    = 100 * mean(min(max(keepRatio(keepZone, :), 0), 1), 1).';
ZeroTap_dB   = 10 * log10(squeeze(mean(zeroEnergy, [1 2]))) ;
ZeroTap_dB   = ZeroTap_dB(:);

%% Display results
fprintf("\n=== Final NMSD (dB) vs small-tap magnitude s (cutoff = %.3g) ===\n\n", cutoff);
disp(array2table(final_dB, 'RowNames', cellstr("s_" + string(st.smallVals)), ...
                 'VariableNames', cellstr(stNames)));

fprintf("\n=== Relative bias of the small tap (%%), mean over channels ===\n");
fprintf("Negative = estimate pulled toward zero\n\n");
disp(array2table(meanBias, 'RowNames', cellstr("s_" + string(st.smallVals)), ...
                 'VariableNames', cellstr(stNames)));

fprintf("\n=== Small tap suppressed (%% of channels, estimate below %.0f%% of s) ===\n\n", ...
        100 * st.suppressFrac);
disp(array2table(suppressPct, 'RowNames', cellstr("s_" + string(st.smallVals)), ...
                 'VariableNames', cellstr(stNames)));

fprintf("\n=== Thresholding view (support recovery) ===\n");
fprintf("Removal score  : %% of a tap removed, averaged over s <= %.3g (should be high)\n", ...
        st.removeBelow * cutoff);
fprintf("Keep score     : %% of a tap kept, averaged over %.3g <= s <= %.3g (should be high)\n", ...
        st.keepFrom * cutoff, cutoff);
fprintf("ZeroTap_dB     : leftover energy on true-zero taps vs ||h||^2 (should be low)\n\n");
threshTable = table(stNames(:), RemovalScore, KeepScore, ZeroTap_dB, ...
    'VariableNames', {'Algorithm', 'Removal_pct', 'Keep_pct', 'ZeroTap_dB'});
disp(threshTable);

%% Plot styles
colors     = turbo(nAlg + 2);
colors     = colors(2:end-1, :);        % skip near-black and near-white ends
lineStyles = ["-", "--", "-.", ":"];

%% Attraction functions (theory)
wGrid = linspace(0, 1.5 * cutoff, 600);
figure;
hold on;
for k = 1:nAlg
    if stNames(k) == "lms_c"
        continue;
    end
    plot(wGrid, attraction(stNames(k), wGrid, alg), ...
        'Color', colors(k, :), ...
        'LineStyle', lineStyles(mod(k - 1, numel(lineStyles)) + 1), ...
        'LineWidth', 1.8);
end
xline(cutoff, 'k:', 'cutoff', 'HandleVisibility', 'off');
grid on;
xlabel('|w|');
ylabel('Attraction c(|w|)');
title('Attraction functions (same cutoff)');
legend(sparse, 'Interpreter', 'none', 'Location', 'northeast');
hold off;

%% Final NMSD vs small-tap magnitude
figure;
hold on;
for k = 1:nAlg
    errorbar(st.smallVals, final_dB(:, k), finalLo_dB(:, k), finalHi_dB(:, k), '-o', ...
        'Color', colors(k, :), ...
        'LineStyle', lineStyles(mod(k - 1, numel(lineStyles)) + 1), ...
        'LineWidth', 1.4);
end
xline(cutoff, 'k:', 'cutoff', 'HandleVisibility', 'off');
grid on;
xlabel('Small tap magnitude s');
ylabel('Final NMSD (dB)');
title('Final NMSD vs. small-tap magnitude (95% CI)');
legend(stNames, 'Interpreter', 'none', 'Location', 'best');
hold off;

%% Small-tap bias: measured vs predicted
figure;
hold on;
for k = 1:nAlg
    errorbar(st.smallVals, meanBias(:, k), ciBias(:, k), '-o', ...
        'Color', colors(k, :), ...
        'LineStyle', lineStyles(mod(k - 1, numel(lineStyles)) + 1), ...
        'LineWidth', 1.4, ...
        'DisplayName', stNames(k) + " (measured)");
end
for k = 1:nAlg
    if stNames(k) == "lms_c"
        continue;
    end
    plot(st.smallVals, predBias(:, k), ':', ...
        'Color', colors(k, :), ...
        'LineWidth', 1.2, ...
        'DisplayName', stNames(k) + " (predicted)");
end
xline(cutoff, 'k:', 'cutoff', 'HandleVisibility', 'off');
yline(0, 'k-', 'HandleVisibility', 'off');
grid on;
xlabel('Small tap magnitude s');
ylabel('Relative bias of small tap (%)');
title('How much each method pulls a real small tap toward zero');
legend('Interpreter', 'none', 'Location', 'southeast');
hold off;

%% Suppression rate
figure;
hold on;
for k = 1:nAlg
    plot(st.smallVals, suppressPct(:, k), '-o', ...
        'Color', colors(k, :), ...
        'LineStyle', lineStyles(mod(k - 1, numel(lineStyles)) + 1), ...
        'LineWidth', 1.4);
end
xline(cutoff, 'k:', 'cutoff', 'HandleVisibility', 'off');
grid on;
xlabel('Small tap magnitude s');
ylabel('Channels with small tap suppressed (%)');
title(sprintf('Small tap estimated below %.0f%% of its true value', 100 * st.suppressFrac));
legend(stNames, 'Interpreter', 'none', 'Location', 'northeast');
hold off;

%% Thresholding view: fraction of the small tap kept vs ideal threshold
figure;
hold on;
yl = [-0.05 1.1];
patch([st.removeBelow st.keepFrom st.keepFrom st.removeBelow] * cutoff, ...
      [yl(1) yl(1) yl(2) yl(2)], [0.5 0.5 0.5], ...
      'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
idealS = [st.smallVals(1), st.removeBelow * cutoff, NaN, st.keepFrom * cutoff, st.smallVals(end)];
idealR = [0, 0, NaN, 1, 1];
plot(idealS, idealR, 'k-', 'LineWidth', 2.5, 'DisplayName', 'ideal threshold');
for k = 1:nAlg
    plot(st.smallVals, keepRatio(:, k), '-o', ...
        'Color', colors(k, :), ...
        'LineStyle', lineStyles(mod(k - 1, numel(lineStyles)) + 1), ...
        'LineWidth', 1.4, ...
        'DisplayName', stNames(k));
end
xline(cutoff, 'k:', 'cutoff', 'HandleVisibility', 'off');
ylim(yl);
grid on;
xlabel('Small tap magnitude s');
ylabel('Fraction of small tap kept, E[w_{small}] / s');
title('Thresholding behaviour (shaded = no preference)');
legend('Interpreter', 'none', 'Location', 'southeast');
hold off;

%% Zero-tap energy vs small-tap magnitude
zeroE_dB = 10 * log10(squeeze(mean(zeroEnergy, 1)));   % nS x nAlg
figure;
hold on;
for k = 1:nAlg
    plot(st.smallVals, zeroE_dB(:, k), '-o', ...
        'Color', colors(k, :), ...
        'LineStyle', lineStyles(mod(k - 1, numel(lineStyles)) + 1), ...
        'LineWidth', 1.4);
end
xline(cutoff, 'k:', 'cutoff', 'HandleVisibility', 'off');
grid on;
xlabel('Small tap magnitude s');
ylabel('Energy on true-zero taps / ||h||^2 (dB)');
title('False-tap energy (lower = cleaner zero taps)');
legend(stNames, 'Interpreter', 'none', 'Location', 'best');
hold off;

%% Save results
save("small_tap_results.mat", "st", "stNames", "finalNMSD", "relBias", ...
     "suppressed", "zeroEnergy", "threshTable", "keepRatio", "final_dB", "meanBias", "predBias", "suppressPct", "alg", "cfg");
fprintf("\nSaved results to small_tap_results.mat\n");

%% Local functions

function c = attraction(name, w, alg)
%ATTRACTION Attraction function c(|w|) of each sparse algorithm (peak 1).

w = abs(w);

switch name
    case "cza_lms"
        a = alg.cza_lms.a;
        c = 0.5 * (1 + cos(pi * w / a)) .* (w <= a);
    case "l0_lms"
        beta = alg.l0_lms.beta;
        c = (1 - beta * w) .* (w <= 1 / beta);
    case "pza_lms"
        alpha = alg.pza_lms.alpha;
        c = (1 - (alpha - 1) * w) ./ (1 + w).^(alpha + 1) .* (w <= 1 / (alpha - 1));
    otherwise
        c = zeros(size(w));
end
end