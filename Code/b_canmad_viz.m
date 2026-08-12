%% ========================================================================
% Project: Exploring CAN MD: Gagnon et al.., Replication
% Author: Amal Varghese
% Purpose: PCA Analysis - Figures 2 and 4
% ========================================================================

%% Prerequisites
clear;
clc;
close all;

%% Load PCA data from a_canmad_pca.m
load(fullfile('Output', 'pca_data.mat'));
disp('PCA data loaded successfully.');
disp(size(pca_matrix_std));

%% FIGURE 2- Eigenvalues and explanatory power of factors
fprintf('Sum of all eigenvalues: %.1f\n', sum(latent));
fprintf('Number of series: %d\n', size(pca_matrix_std, 2));

n_comp      = 20;
eigenvalues = latent(1:n_comp);
cumvar      = cumsum(explained(1:n_comp));
fprintf('First eigenvalue: %.2f\n', eigenvalues(1));
fprintf('Cumulative variance at PC10: %.1f%%\n', cumvar(10));

% Figure 2
fig2 = figure;
t = tiledlayout(2,1,'TileSpacing', 'compact','Padding','compact');
nexttile;  % Fig 2 panel a
plot(1:n_comp,eigenvalues, 'k-','LineWidth', 1.5);
ylabel('Eigenvalues');
xlabel('Number of factors');
title('Scree plot');
nexttile; % Fig 2 panel b
plot(1:n_comp, cumvar, 'k-', 'LineWidth', 1.5);
ylabel('% of variance explained');
xlabel('Number of factors');
title('Trace');
nexttile(1); % modifying X axis
xlim([1 n_comp]);
nexttile(2);
xlim([1 n_comp]);
nexttile(1); %adding grid lines
grid on;
nexttile(2);
grid on;


%% FIGURE 4- Factors 1 to 4 and their main series
% Identify the highest loading variable for each of the first 4 factors
[~, top_idx] = max(abs(coeff(:, 1:4)));
top_series = aggregate_vars(top_idx);
disp(top_series');
% dates
dates = pca_df.Date(1:end-1);
fprintf('Date range: %s to %s\n', datestr(dates(1)), datestr(dates(end)));
% figure
fig4 = figure;
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile; % factor 1
hold on;
patch([datetime('2008-10-01') datetime('2009-05-01') datetime('2009-05-01') datetime('2008-10-01')], ...
      [-6 -6 6 6], [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
yyaxis left;
plot(dates, score(:,1) ./ std(score(:,1)), 'k-', 'LineWidth', 1);
ylabel('First factor');
yyaxis right;
plot(dates, pca_df{1:end-1, top_series{1}} ./ std(pca_df{1:end-1, top_series{1}}), 'r-', 'LineWidth', 1);
ylabel(top_series{1});
title('Factor 1, GDP');

nexttile; % factor 2
hold on;
patch([datetime('2008-10-01') datetime('2009-05-01') datetime('2009-05-01') datetime('2008-10-01')], ...
      [-6 -6 6 6], [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
yyaxis left;
plot(dates, score(:,2) ./ std(score(:,2)), 'k-', 'LineWidth', 1);
ylabel('Second factor');
yyaxis right;
plot(dates, pca_df{1:end-1, top_series{2}} ./ std(pca_df{1:end-1, top_series{2}}), 'r-', 'LineWidth', 1);
ylabel(top_series{2});
title('Factor 2, Exchange Rate');

nexttile; % factor 3
hold on;
patch([datetime('2008-10-01') datetime('2009-05-01') datetime('2009-05-01') datetime('2008-10-01')], ...
      [-6 -6 6 6], [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
yyaxis left;
plot(dates, score(:,3) ./ std(score(:,3)), 'k-', 'LineWidth', 1);
ylabel('Third factor');
yyaxis right;
plot(dates, pca_df{1:end-1, top_series{3}} ./ std(pca_df{1:end-1, top_series{3}}), 'r-', 'LineWidth', 1);
ylabel(top_series{3});
title('Factor 3, Commercial Paper Rate');

nexttile; % factor 4
hold on;
patch([datetime('2008-10-01') datetime('2009-05-01') datetime('2009-05-01') datetime('2008-10-01')], ...
      [-6 -6 6 6], [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.4);
yyaxis left;
plot(dates, score(:,4) ./ std(score(:,4)), 'k-', 'LineWidth', 1);
ylabel('Fourth factor');
yyaxis right;
plot(dates, pca_df{1:end-1, top_series{4}} ./ std(pca_df{1:end-1, top_series{4}}), 'r-', 'LineWidth', 1);
ylabel(top_series{4});
title('Factor 4, Bond Spread (5-10Y)');

% additional cosmetics
% axis color
nexttile(1); ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';
nexttile(2); ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';
nexttile(3); ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';
nexttile(4); ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'k';

% gridlines
nexttile(1); grid on;
nexttile(2); grid on;
nexttile(3); grid on;
nexttile(4); grid on;

% Save Plots
exportgraphics(figure(1), fullfile('Output', 'Figure2_scree.pdf'), 'ContentType', 'vector');
exportgraphics(figure(2), fullfile('Output', 'Figure4_factors.pdf'), 'ContentType', 'vector');