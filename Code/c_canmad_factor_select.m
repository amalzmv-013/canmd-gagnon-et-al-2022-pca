%% ========================================================================
% Project: Exploring CAN MD: Gagnon et al.., Replication
% Author: Amal Varghese
% Purpose: PCA Analysis - Table 1 and Figure 3
% ========================================================================

%% Prerequisites
clear;
clc;
close all;

%% Load PCA data from a_canmad_pca.m
load(fullfile('Output', 'pca_data.mat'));
disp('PCA data loaded successfully.');
fprintf('Matrix dimensions: %d rows x %d columns\n', size(pca_matrix_std,1), size(pca_matrix_std,2));

%% Adding necessary tool box made by the author to the path
addpath(genpath(fullfile('Tools','EstimationNumberFactors')));
disp('Toolbox loaded successfully.');

%% Load balanced CAN-MD panel [EM-imputed version with no missing values]
data_can = readtable (fullfile('Input','balanced_can_md.csv'));
fprintf('Balanced panel dimensions: %d rows x %d columns\n', height(data_can), width(data_can));
bal_vars = data_can.Properties.VariableNames;
fprintf('Total columns including Date: %d\n', numel(bal_vars));

%% Filter to aggregate series only
bal_vars_noDate = setdiff(bal_vars, {'Date'}, 'stable');
prov_pattern = '_(NF|PEI|NS|NB|QC|ONT|MAN|SAS|ALB|BC)(_|$)';
is_prov = ~cellfun(@isempty, regexp(bal_vars_noDate, prov_pattern, 'once'));
agg_vars_bal = bal_vars_noDate(~is_prov);
fprintf('Aggregate series retained: %d\n', numel(agg_vars_bal));

%% Extract numeric matrix
data_can_agg = table2array(data_can(:, agg_vars_bal));
fprintf('Matrix dimensions: %d rows x %d columns\n', size(data_can_agg, 1), size(data_can_agg, 2));
fprintf('NaNs remaining: %d\n', sum(isnan(data_can_agg(:))));
%% Check select_numfac input expectations
help select_numfac; %[two inputs are to be provided X and rmax]
open select_numfac;

%% TABLE 1 Estimate number of factors 
rmax = 20;
fprintf('Running select_numfac... this may take a moment.\n');
estim_numfac = select_numfac(data_can_agg, rmax);
disp('Done.');
% Display  estimated factors
BN02 = estim_numfac.bn2002(2,2);
ABC = estim_numfac.abc;
ON  = estim_numfac.o2010;
AH  = estim_numfac.ah2013(1,1);
HL  = estim_numfac.hl2007;
BN07 = estim_numfac.b2007(1,1);
AW  = estim_numfac.aw2007(2,2);
fprintf('BN02: %d factors\n', BN02);
fprintf('ABC:  %d factors\n', ABC);
fprintf('ON:   %d factors\n', ON);
fprintf('AH:   %d factors\n', AH);
fprintf('HL:   %d factors\n', HL);
fprintf('BN07: %d factors\n', BN07);
fprintf('AW:   %d factors\n', AW);
% Display Table 1
criteria = {'BN02'; 'ABC'; 'ON'; 'AH'; 'HL'; 'BN07'; 'AW'};
estimates = [BN02; ABC; ON; AH; HL; BN07; AW];
Table1 = table(criteria, estimates, 'VariableNames', {'Criterion', 'N_Factors'});
disp('Table 1 — Number of Factors (Canada panel):');
disp(Table1);
% save
writetable(Table1, fullfile('Output', 'Table1_factors.csv'));
fprintf('Table 1 saved.\n');

%% Figure 3 — Recursive factor selection
% Define recursive window
t1 = (1991 - 1981) * 12 + 1;
TT = height(data_can_agg) - t1 + 1;

% Pre-allocate results
BN2002_rec = zeros(TT, 1);
HL2007_rec = zeros(TT, 1);

% Recursive expanding window loop
% Note: HL must run before BN02 within each iteration due to legacy
% random number generator conflict inside HLestimate()
fprintf('Running recursive factor selection... this will take a few minutes.\n');
for t = 1:TT
    XX = standard(data_can_agg(1:t1+t-1, :));
    HL2007_rec(t) = HLestimate(XX', rmax);
    BN2002_rec(t) = nbplog(XX, rmax, 2, 0);
end
fprintf('Complete. %d iterations.\n', TT);

%% Build date axis
rec_dates = data_can.Date(t1:end);

% Plot Figure 3
fig3 = figure;
tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
hold on;
h1 = plot(rec_dates, BN2002_rec, 'k-',  'LineWidth', 1.5);
h2 = plot(rec_dates, HL2007_rec, 'b--', 'LineWidth', 1.5);
legend([h1 h2], {'Bai-Ng 2002', 'Hallin-Liska 2007'}, 'Location', 'northwest');
ylabel('Number of factors');
xlabel('Time');
ylim([2 6]);
yticks(2:0.5:6);
xlim([rec_dates(1) rec_dates(end)]);
xticks(datetime(1993,1,1):calyears(4):datetime(2017,1,1));
xtickformat('yyyy');
grid on;

%% Save
exportgraphics(fig3, fullfile('Output', 'Figure3_recursive.pdf'), 'ContentType', 'vector');
fprintf('Figure 3 saved.\n');