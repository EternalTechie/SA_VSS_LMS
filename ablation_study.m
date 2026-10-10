
clear;
clc;
close all;

%% Initialize the main simulation
run('main.m');

%% Cutoff parameter values to test
aValues = [0.01 0.03 0.05 0.1 0.15 0.2 0.3 0.5 0.8];
numA = length(aValues);

%% Storage
czaResults = cell(numA, 1);
czaMetrics = cell(numA, 1);

finalNMSD_dB = zeros(numA, 1);
finalMSD_dB  = zeros(numA, 1);
finalMSE_dB  = zeros(numA, 1);

%% Run CZA-LMS for each cutoff parameter
for k = 1:numA

    algCZA = alg.cza_lms;
    algCZA.a = aValues(k);

    fprintf('Running CZA-LMS with a = %.4f\n', aValues(k));

    czaResults{k} = cza_lms(data, cfg, algCZA);
    czaMetrics{k} = evaluateResult(data, czaResults{k}, cfg);

    finalNMSD_dB(k) = czaMetrics{k}.finalNMSD_dB;
    finalMSD_dB(k)  = czaMetrics{k}.finalMSD_dB;
    finalMSE_dB(k)  = czaMetrics{k}.finalMSE_dB;
end

%% Plot NMSD convergence for all cutoff values
figure;
hold on;
grid on;

for k = 1:numA
    plot(czaMetrics{k}.nmsd_dB, ...
        'DisplayName', sprintf('a = %.3g', aValues(k)));
end

xlabel('Iteration');
ylabel('NMSD (dB)');
title('CZA-LMS: Effect of Cutoff Parameter a');
legend('show', 'Location', 'best');
hold off;

%% Plot final NMSD against a
figure;
semilogx(aValues, finalNMSD_dB, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Cutoff parameter a');
ylabel('Final NMSD (dB)');
title('CZA-LMS: Final NMSD vs. a');

%% Plot final MSD against a
figure;
semilogx(aValues, finalMSD_dB, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Cutoff parameter a');
ylabel('Final MSD (dB)');
title('CZA-LMS: Final MSD vs. a');

%% Plot final MSE against a
figure;
semilogx(aValues, finalMSE_dB, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Cutoff parameter a');
ylabel('Final MSE (dB)');
title('CZA-LMS: Final MSE vs. a');

%% Display results
resultsTable = table( ...
    aValues(:), ...
    finalNMSD_dB, ...
    finalMSD_dB, ...
    finalMSE_dB, ...
    'VariableNames', ...
    {'a', 'FinalNMSD_dB', 'FinalMSD_dB', 'FinalMSE_dB'});

disp(resultsTable);

%% Find the best tested a based on final NMSD
[~, bestIdx] = min(finalNMSD_dB);

fprintf('\nBest tested cutoff parameter a = %.4f\n', ...
    aValues(bestIdx));
fprintf('Final NMSD = %.4f dB\n', finalNMSD_dB(bestIdx));
