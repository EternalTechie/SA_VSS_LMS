%% Rho Sweep
% Sweeps the zero-attraction strength rho for each sparse algorithm over
% many random sparse channels with random sparsity, and saves the best
% rho per algorithm for monte_carlo_study.m.
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

%% Algorithms and rho values to test
sparseNames  = ["za_lms", "rza_lms", "l0_lms", "pza_lms", "cza_lms"];
baselineName = "lms_c";                       % same mu, no attraction
rhoValues    = [1e-5 3e-5 1e-4 3e-4 1e-3];

%% Monte Carlo parameters
rs.N          = 20000;      % iterations per channel
rs.Nch        = 200;        % number of random channels
rs.Kmin       = 2;          % minimum number of active taps
rs.maxDensity = 0.25;       % maximum fraction of active taps
rs.hMin       = 0.2;        % minimum active tap magnitude (threshold)
rs.hMax       = 1.0;        % maximum active tap magnitude
rs.tailFrac   = 0.1;        % final NMSD = mean over last 10% of iterations
rs.seedBase   = 5e5;        % different channels from monte_carlo_study.m

% Convergence-speed metric (computed on the averaged NMSD curve)
rs.target_dB = -30;         % T_target: first iteration below -30 dB
rs.smoothWin = 50;          % moving-average window before thresholding

%% Apply simulation settings
cfg.N = rs.N;
cfg.channel.normalize = false;    % NMSD is already normalized by ||h||^2
cfg.channel.factors   = {};       % static channel

rs.Kmax = max(rs.Kmin, floor(rs.maxDensity * cfg.M));

%% Storage
nAlg = numel(sparseNames);
nRho = numel(rhoValues);

finalNMSD = nan(rs.Nch, nRho, nAlg);    % per-channel final NMSD (linear)
sumNMSD   = zeros(rs.N, nRho, nAlg);    % sum of linear NMSD curves
baseFinal = nan(rs.Nch, 1);
baseSum   = zeros(rs.N, 1);
Kch       = nan(rs.Nch, 1);

tailIdx = floor((1 - rs.tailFrac) * rs.N) + 1 : rs.N;

%% Run sweep
tStart = tic;

for ch = 1:rs.Nch

    % Same channel, input, and noise for every algorithm and rho
    rng(rs.seedBase + ch);

    [taps, gains, K] = drawSparseChannel(cfg.M, rs);
    cfg.channel.activeTaps = taps;
    cfg.channel.gains      = gains;
    Kch(ch) = K;

    d    = cfg.signal.generator(cfg.N, cfg.signal);
    data = generateReceivedData(d, cfg);
    hPow = sum(data.h.^2, 1);

    % Baseline without zero attraction
    res  = runEqualizer(baselineName, data, cfg, alg);
    nmsd = (sum((res.w - data.h).^2, 1) ./ hPow).';
    baseSum       = baseSum + nmsd;
    baseFinal(ch) = mean(nmsd(tailIdx));

    % Sparse algorithms at every rho
    for k = 1:nAlg
        algK = alg;
        for j = 1:nRho
            algK.(sparseNames(k)).rho = rhoValues(j);

            res  = runEqualizer(sparseNames(k), data, cfg, algK);
            nmsd = (sum((res.w - data.h).^2, 1) ./ hPow).';

            sumNMSD(:, j, k)    = sumNMSD(:, j, k) + nmsd;
            finalNMSD(ch, j, k) = mean(nmsd(tailIdx));
        end
    end

    if mod(ch, 10) == 0
        fprintf("Channel %d/%d | %.0f s elapsed\n", ch, rs.Nch, toc(tStart));
    end

end

%% Results
meanFinal_dB = reshape(10 * log10(mean(finalNMSD, 1)), nRho, nAlg);
baseFinal_dB = 10 * log10(mean(baseFinal));

T_target = nan(nRho, nAlg);
for k = 1:nAlg
    for j = 1:nRho
        T_target(j, k) = targetTime(sumNMSD(:, j, k) / rs.Nch, rs);
    end
end
baseT_target = targetTime(baseSum / rs.Nch, rs);

[bestFinal_dB, bestIdx] = min(meanFinal_dB, [], 1);

bestRho = struct();
for k = 1:nAlg
    bestRho.(sparseNames(k)) = rhoValues(bestIdx(k));
    if bestIdx(k) == 1 || bestIdx(k) == nRho
        warning("Best rho for %s is at the edge of the grid (%.3g). " + ...
                "Extend rhoValues.", sparseNames(k), rhoValues(bestIdx(k)));
    end
end

%% Display results
fprintf("\n=== Final NMSD (dB) vs rho, %d channels, K = %d to %d ===\n", ...
        rs.Nch, rs.Kmin, rs.Kmax);
fprintf("Baseline %s: %.3f dB, T_target = %d\n\n", ...
        baselineName, baseFinal_dB, baseT_target);

finalTable = array2table(meanFinal_dB, ...
    'RowNames', cellstr("rho_" + string(rhoValues)), ...
    'VariableNames', cellstr(sparseNames));
disp(finalTable);

fprintf("\n=== Iterations to reach %g dB vs rho ===\n\n", rs.target_dB);
targetTable = array2table(T_target, ...
    'RowNames', cellstr("rho_" + string(rhoValues)), ...
    'VariableNames', cellstr(sparseNames));
disp(targetTable);

fprintf("\n=== Best rho per algorithm ===\n\n");
for k = 1:nAlg
    j = bestIdx(k);
    fprintf("%-8s | best rho = %.3g | final NMSD = %.3f dB | gain over %s = %.2f dB | T_target = %d\n", ...
            sparseNames(k), rhoValues(j), bestFinal_dB(k), baselineName, ...
            baseFinal_dB - bestFinal_dB(k), T_target(j, k));
end

%% Plot styles
colors     = turbo(nAlg + 2);
colors     = colors(2:end-1, :);        % skip near-black and near-white ends
lineStyles = ["-", "--", "-.", ":"];
plotIdx    = unique([1:10:rs.N, rs.N]);

%% Final NMSD vs rho
figure;
hold on;
for k = 1:nAlg
    plot(rhoValues, meanFinal_dB(:, k), '-o', ...
        'Color', colors(k, :), ...
        'LineStyle', lineStyles(mod(k - 1, numel(lineStyles)) + 1), ...
        'LineWidth', 1.5);
end
yline(baseFinal_dB, 'k:', baselineName, 'Interpreter', 'none');
set(gca, 'XScale', 'log');
grid on;
xlabel('\rho');
ylabel('Final NMSD (dB)');
title(sprintf('Final NMSD vs. \\rho (%d channels)', rs.Nch));
legend([sparseNames, baselineName], 'Interpreter', 'none', 'Location', 'best');
hold off;

%% Convergence speed vs rho
figure;
hold on;
for k = 1:nAlg
    plot(rhoValues, T_target(:, k), '-o', ...
        'Color', colors(k, :), ...
        'LineStyle', lineStyles(mod(k - 1, numel(lineStyles)) + 1), ...
        'LineWidth', 1.5);
end
yline(baseT_target, 'k:', baselineName, 'Interpreter', 'none');
set(gca, 'XScale', 'log');
grid on;
xlabel('\rho');
ylabel(sprintf('Iterations to reach %g dB', rs.target_dB));
title('Convergence speed vs. \rho');
legend([sparseNames, baselineName], 'Interpreter', 'none', 'Location', 'best');
hold off;

%% Averaged NMSD at best rho
figure;
hold on;
for k = 1:nAlg
    j = bestIdx(k);
    plot(plotIdx, 10 * log10(sumNMSD(plotIdx, j, k) / rs.Nch), ...
        'Color', colors(k, :), ...
        'LineStyle', lineStyles(mod(k - 1, numel(lineStyles)) + 1), ...
        'LineWidth', 1.4);
end
plot(plotIdx, 10 * log10(baseSum(plotIdx) / rs.Nch), 'k:', 'LineWidth', 1.4);
grid on;
xlabel('Iteration');
ylabel('Average NMSD (dB)');
title('Average NMSD at best \rho per algorithm');
legend([sparseNames + " (rho=" + string(rhoValues(bestIdx)) + ")", baselineName], ...
       'Interpreter', 'none', 'Location', 'best');
hold off;

%% Save results
save("rho_sweep_results.mat", "bestRho", "rhoValues", "sparseNames", ...
     "meanFinal_dB", "T_target", "baseFinal_dB", "baseT_target", "Kch", "rs", ...
     "alg", "cfg");
fprintf("\nSaved best rho values to rho_sweep_results.mat\n");

%% Local functions

function [taps, gains, K] = drawSparseChannel(M, rs)
%DRAWSPARSECHANNEL Random sparse channel: K active taps at random positions,
% each with magnitude drawn directly from [rs.hMin, rs.hMax] and random sign.

K     = randi([rs.Kmin, rs.Kmax]);
taps  = sort(randperm(M, K));
mags  = rs.hMin + (rs.hMax - rs.hMin) * rand(K, 1);
signs = 2 * (rand(K, 1) >= 0.5) - 1;
gains = mags .* signs;
end

function tTarget = targetTime(avgLin, rs)
%TARGETTIME First iteration where the smoothed averaged NMSD reaches the
% target (NaN if never reached).

smooth_dB = movmean(10 * log10(avgLin), rs.smoothWin);
tTarget = find(smooth_dB <= rs.target_dB, 1);
if isempty(tTarget)
    tTarget = NaN;
end
end
