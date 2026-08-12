%% ========================================================================
% Project: Exploring CAN MD: Gagnon et al.., Replication
% Author: Amal Varghese
% Purpose: PCA Analysis (Preliminary Set up) 
% ========================================================================

%% Prerequisites
clear;
clc;
close all;

%% Load and transfrom the CANMD database
file_path = fullfile('Input', 'TR_CAN_MD.csv');
can_md_raw = readtable(file_path);

% Defining boundary dates
startDate = datetime('1992-01-01');
endDate = datetime('2020-01-01');
can_md_filtered = can_md_raw(can_md_raw.Date >= startDate & can_md_raw.Date <= endDate, :);

% Summary Stat
key_vars = can_md_filtered(:, {'Date', 'GDP_new', 'EMP_CAN', 'Imp_BP_new', 'Exp_BP_new', 'CPI_ALL_CAN', 'SP500'});
key_vars = convertvars(key_vars, {'GDP_new', 'EMP_CAN', 'Imp_BP_new', 'Exp_BP_new', 'SP500'}, @str2double);
summary(key_vars);

% Correlation on Key variables
data_matrix = table2array(key_vars(:, 2:end));
corr_matrix = corrcoef(data_matrix, 'Rows', 'complete');
disp(corr_matrix);

% Plots
% Create a 2*3 panel plot of the key_vars
figure;
tiledlayout(3,2);

nexttile;
plot(key_vars.Date, key_vars.GDP_new);
title('GDP');

nexttile;
plot(key_vars.Date, key_vars.EMP_CAN);
title('Employment');

nexttile;
plot(key_vars.Date, key_vars.Imp_BP_new);
title('Imports (BoP)');

nexttile;
plot(key_vars.Date, key_vars.Exp_BP_new);
title('Exports (BoP)');

nexttile;
plot(key_vars.Date, key_vars.CPI_ALL_CAN);
title('CPI (All-items)');

nexttile;
plot(key_vars.Date, key_vars.SP500);
title('S&P 500');

%% PCA DataSet Creation

% Screening the column names of the data set
all_vars = can_md_filtered.Properties.VariableNames;
all_vars_noDate = setdiff(all_vars, {'Date'}, 'stable');
disp(all_vars(1:5));
fprintf('Total series (excl. Date): %d\n', numel(all_vars_noDate));

% Screening for column names with provincial suffix
prov_pattern = '_(NF|PEI|NS|NB|QC|ONT|MAN|SAS|ALB|BC)(_|$)';
disp(prov_pattern);
match_result = regexp(all_vars_noDate, prov_pattern, 'once');
disp(match_result(1:5))
is_provincial = ~cellfun(@isempty, match_result);
fprintf('Provincial series flagged: %d\n', sum(is_provincial));

% Filtering out provincial series columns
aggregate_vars = all_vars_noDate(~is_provincial);
fprintf('Aggregate series retained: %d\n', numel(aggregate_vars));

% Slicing CAN-MD for aggregate variables from 1992-2020
pca_df = can_md_filtered(:, ['Date', aggregate_vars]);
fprintf('pca_df dimensions: %d rows x %d columns\n', height(pca_df), width(pca_df));

%Converting cell-type columns 
var_types = varfun(@class, pca_df(:, 2:end), 'OutputFormat', 'cell');
is_cell_var = strcmp(var_types, 'cell');
cell_vars = aggregate_vars(is_cell_var);
fprintf('Columns needing conversion: %d\n', numel(cell_vars));
pca_df = convertvars(pca_df, cell_vars, @str2double);
pca_matrix = table2array(pca_df(:, 2:end));
fprintf('pca_matrix dimensions: %d rows x %d columns\n', size(pca_matrix, 1), size(pca_matrix, 2));

% Adressing NaN
nan_counts = sum(isnan(pca_matrix), 1);
fprintf('Columns with any NaN: %d\n', sum(nan_counts > 0));
fprintf('Columns fully clean: %d\n', sum(nan_counts == 0));
nan_series = aggregate_vars(nan_counts > 0);
disp(nan_series');
nan_counts_flagged = nan_counts(nan_counts > 0);
disp(table(nan_series', nan_counts_flagged', 'VariableNames', {'Series', 'NaN_Count'}));
[nan_rows, nan_cols] = find(isnan(pca_matrix));
disp(table(nan_rows, aggregate_vars(nan_cols)', 'VariableNames', {'RowIndex', 'Series'}));
pca_matrix = pca_matrix(:, ~strcmp(aggregate_vars, 'FIN_new'));
aggregate_vars = aggregate_vars(~strcmp(aggregate_vars, 'FIN_new'));
fprintf('Dimensions after dropping FIN_new: %d rows x %d columns\n', size(pca_matrix,1), size(pca_matrix,2));
pca_matrix = pca_matrix(1:end-1, :);
fprintf('Dimensions after dropping last row: %d rows x %d columns\n', size(pca_matrix,1), size(pca_matrix,2));
fprintf('Remaining NaNs: %d\n', sum(isnan(pca_matrix(:))));

% Standardize the matrix
pca_matrix_std = zscore(pca_matrix);
fprintf('Column 1 mean: %.6f\n', mean(pca_matrix_std(:,1)));
fprintf('Column 1 std:  %.6f\n', std(pca_matrix_std(:,1)));

%% PCA 
[coeff, score, latent, ~, explained] = pca(pca_matrix_std, 'Centered', false);
fprintf('Number of components extracted: %d\n', size(coeff, 2));
for i = 1:10
    fprintf('PC%d: %.1f%% variance explained (Cumulative: %.1f%%)\n', ...
        i, explained(i), sum(explained(1:i)));
end 

%% Save
save(fullfile('Output', 'pca_data.mat'), 'pca_matrix_std', 'coeff', 'score', 'latent', ...
     'explained', 'aggregate_vars', 'pca_df', 'corr_matrix');