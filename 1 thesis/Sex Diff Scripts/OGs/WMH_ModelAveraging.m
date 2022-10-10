%% ************************* Model Averaging ****************************
%
% Lisa Röhrig, 06-29-2021, last update on 09-08-2021
%
% For final model performance, average predictions of multiple model 
% repetitions that differ due to randomization of CV-folds
%
% >> Models:
%       * Model 1 uses Lesion data 
%       * Model 2 uses Lesion + WMH data
%       * Model 3 uses Lesion + permuted WMH data
%       * Model 4 uses Lesion data (+ WMH data) + age
%
% >> WMH data formats:
%       * Type 1 uses volumetric data of the whole brain
%       * Type 2 uses volumetric data of the left hemisphere
%       * Type 3 uses volumetric data of the right hemisphere
%       * Type 4 uses rating scores (CHS/Fazekas)

clear; clc; close all;

% Path requirements for this script: within the main folder, there should be folders for
% structural data ('Data') and for output ('Results'). The excel-file with WMH
% ratings, behavior etc. should be within the 'Data' folder. The results-folder 
% should contain folders for each model version (e.g., 'Model_1', 'Model_2\Addition', 
% 'Model_3\Concatenation'), each of those subfolders should contain an additional 
% 'ModelAveraging'-folder.  
path = 'D:\Lisa\WMH'; % Path needs to be modified
cd(path);

%% Definition of parameters 
fprintf('***Model Averaging ***\n');

% Define mode
fprintf('Model versions: 1 - Lesion data, 2 - Lesion & WMH data, 3 - Lesion & permuted WMH data, 4 - Lesion (& WMH data) & age\n');
model_version = input('Select model version. Version # ');
while (model_version ~= 1) && (model_version ~= 2) && (model_version ~= 3) && (model_version ~= 4)
   model_version = input('Select model version. Number should be between 1 and 4! Version # ');
end
if model_version == 4 % Age
    age_mod = input('Age model (1) with or (2) without WMH data? ...');
    while (age_mod ~= 1) && (age_mod ~= 2)
        age_mod = input('Age model (1) with or (2) without WMH data? Number should be 1 or 2! ...');
    end
else 
    age_mod = 0;
end

% Define data type
voxels = 3; % Default (no WMH voxels)
if (model_version ~= 1) && (age_mod ~= 2) % Only if WMH data
    fprintf('Data formats: \t*Volumetric data*: 1 - Whole brain, 2 - Left HS, 3 - Right HS; \t*Rating scores*: 4 - Whole brain\n');
    data_format = input('Select data format. Format # ');
    while (data_format ~= 1) && (data_format ~= 2) && (data_format ~= 3) && (data_format ~= 4)
       data_format = input('Select data format. Number should be between 1 and 4! Format # ');
    end
    % Define whether to use all WMH voxels or voxels that are damaged in at least 
    % five patients
    if data_format ~= 4 % Not rating scores
        voxels = input('Select variant of WMH voxel selection: 1 - all voxels, 2 - ignoring rarely damaged voxels. #');
            while (voxels ~= 1) && (voxels ~= 2)
                voxels = input('Select variant of WMH voxel selection. Number should be 1 or 2! # ');
            end
    else
        % Define which rating scale to use
        fprintf('Rating scales:\t1 - CHS, 2 - Fazekas PV-WMH, 3 - Fazekas DS-WMH, 4 - PV-WMH & DS-WMH');
        rating = input('Select rating scale. Rating # ');
        while (rating ~= 1) && (rating ~= 2) && (rating ~= 3) && (rating ~= 4)
            rating = input('Select rating scale. Number should be 1 or 2! # ');
        end
    end
elseif (model_version == 1) || (model_version == 4)
    data_format = 1;
end

% Define fusion type of lesion + WMH matrices
if (model_version ~= 1)  && (data_format ~= 4) && (age_mod ~= 2)
    mod_var = input('Define type of matrices fusion (1 - matrix concatenation, 2 - matrix addition). ... ');
    while (mod_var ~= 1) && (mod_var ~= 2)
            mod_var = input('Define type of matrices fusion (1 - matrix concatenation, 2 - matrix addition). ... ');
    end
    switch mod_var
        case 1
            mod_txt = 'Conc_';
            mod_folder = 'Concatenation';
        case 2
            mod_txt = 'Add_';
            mod_folder = 'Addition';
    end
elseif model_version == 1 % Lesion
    mod_txt = []; mod_folder = [];
elseif data_format == 4
    mod_txt = 'Conc_'; 
    mod_folder = 'Concatenation';
end
if model_version == 4
    mod_folder = [];
end
switch data_format % For structure name later
    case 1
        if age_mod == 2
            data_txt = 'Conc_Age';
        else
            data_txt = [mod_txt 'VolWB'];
        end
    case 2
        data_txt = [mod_txt 'VolLHS'];
    case 3
        data_txt = [mod_txt 'VolRHS'];
    case 4
        if rating == 1
            data_txt = [mod_txt 'CHS'];
        elseif rating == 2
            data_txt = [mod_txt 'PVWMH'];
        elseif rating == 3
            data_txt = [mod_txt 'DSWMH'];
        elseif rating == 4
            data_txt = [mod_txt 'Fazekas'];
        end
end
if (model_version == 4) && (age_mod ~= 2)
    data_txt = [data_txt '_Age'];
end
switch voxels 
    case 1
        vox_txt = 'AllVox';
    case 2
        vox_txt = 'ShortVox';
    case 3
        vox_txt = '';
end
vox_struc{1} = 'AllVox';
vox_struc{2} = 'ShortVox';
vox_struc{3} = '';

% Number of permutations
if model_version == 3
    n_perm = input('Please enter number of applied permutation-iterations (default = 5000). n = ...\n');
    if isempty(n_perm) == 1
        n_perm = 5000; % Default
    end
end

% Number of models to average
num_mod = input('Define number of models you want to use for averaging (Default = 10). ... ');
if isempty(num_mod) == 1 || (abs(round(num_mod)-num_mod) == 0) == 0 
    num_mod = 10; % Default
end

%% Averaging
% Model path
model_path = fullfile(path, 'Results', ['Model_' num2str(model_version)], mod_folder);
d = dir(model_path);

% Load test outcome for permutation averaging
% Select excel-file containing data 
if model_version == 3
    [num, txt, raw] = xlsread(fullfile(path, 'Data', 'WMH_Data'), 'WMH_Final');
    if isequal(num(:, 1), sort(num(:, 1))) == 0
        fprintf('\nError! Data file has to be sorted regarding RHLM-number!!!\n\n');
        return
    end
    outcome_org = num(:, 20); % Mean CoC 
    outcome_pos = outcome_org; 
    outcome_pos(outcome_pos < 0) = 0; % Only positive CoC (negative ones -> 0)
    test_sort = realsqrt(outcome_pos);
end

% Load structures
fprintf('\nLoading models ....\n');
for n = 1:num_mod
    for i = 1:length(d)
        if model_version ~= 3
            if (contains(d(i).name, 'model_struc_') == 1) && (contains(d(i).name, data_txt) == 1) ...
                    && (contains(d(i).name, ['_rs_' num2str(n) '_']) == 1) ...
                    && (contains(d(i).name, vox_struc{voxels}) == 1) % Select model file with specific name
                curr_mod = d(i).name;
                curr_int = extractAfter(curr_mod, '_avRSquare_');
                RSquare_curr = extractBefore(curr_int, length(curr_int)-21);
                fprintf('Random stream number %d: RSquare = %s\n', n, RSquare_curr);
                break
            end
        elseif model_version == 3
            if (contains(d(i).name, ['Permutation_predictions_' num2str(n_perm)]) == 1) && (contains(d(i).name, data_txt) == 1) ...
                    && (contains(d(i).name, ['_rs_' num2str(n) '_']) == 1) ...
                    && (contains(d(i).name, vox_struc{voxels}) == 1)
                curr_mod = d(i).name;
                curr_int = extractAfter(curr_mod, '_Fit_');
                RSquare_curr = extractBefore(curr_int, length(curr_int)-26);
                fprintf('Random stream number %d: RSquare (0.95-best model) = %s\n', n, RSquare_curr);
                break
            end
        end
    end
    load([model_path '\' curr_mod]);
    
    % Load predictions of nested CVs
    if model_version ~= 3
        pred_all = []; test_all = []; subj_all = [];
        for s = 1:size(model_struc, 2)
            pred_all = [pred_all; model_struc(s).predicted_label_oL];
            test_all = [test_all; model_struc(s).outcome_test_oL];
            subj_all = [subj_all; model_struc(s).subjects_test_oL];
        end
    
        % Sort order of subjects (and predictions)
        [subj_sort, sort_ind] = sort(subj_all);
        pred_sort = pred_all(sort_ind, :);
        test_sort = test_all(sort_ind, :);

        % Save predictions for all models
        mat_all(n).pred_all = pred_all;
        mat_all(n).test_all = test_all;
        mat_all(n).subj_all = subj_all;
        mat_sort(n).pred_sort = pred_sort;
        mat_sort(n).test_sort = test_sort;
        mat_sort(n).subj_sort = subj_sort;

        % Check for same test outcome order
        if n ~= 1
           if isequal(mat_sort(n).test_sort, mat_sort(n-1).test_sort) == 0
               fprintf('\nError! Test outcome must be equal across models to average!\n\n');
               return
           end
        end
    
    % Save sorted predictions of permutations
    elseif model_version == 3
        mat_perm(:, :, n) = perm_pred;
    end
end

% Calculate mean values of predictions
fprintf('\nAveraging ....\n');
if model_version ~= 3
   pred_mean = mean([mat_sort.pred_sort], 2);
elseif model_version == 3
   perm_mean = mean(mat_perm, 3); % calculate mean of all matrix elements across random streams 
end

% Calculate RSquare average
if model_version ~= 3
    squRes = (test_sort - pred_mean).^2;
    squTotal = (test_sort - mean(test_sort)).^2;
    RSquare_av = 1 - (sum(squRes)/sum(squTotal));
    fprintf('\nRSquare-Average (%d random streams): %f\n', num_mod, RSquare_av);
elseif model_version == 3
    for t = 1:size(perm_mean, 2)
        squRes = (test_sort - perm_mean(:, t)).^2;
        squTotal = (test_sort - mean(test_sort)).^2;
        RSquare_mat(t, 1) = 1 - (sum(squRes)/sum(squTotal));
    end
    histogram(RSquare_mat)
    % Sort list 
    perm_sort_fit = sort(RSquare_mat, 1, 'descend');
    % Selecting 95%-best model
    num_p = length(perm_sort_fit) - (length(perm_sort_fit) * 0.95);
    perm_selected = perm_sort_fit(ceil(num_p), :);
    fprintf('\nRSquare of 0.95-best model of averaged permutation predictions (%d random streams): %f\n', num_mod, perm_selected);
end

% Calculate p-value
if model_version == 3    
    % Get model fit of original model
    tmp = dir(fullfile(path, 'Results', ['Model_' num2str(2)], mod_folder, 'ModelAveraging'));
    for m = 1:length(tmp)
        if (contains(tmp(m).name, ['ModelAveraging_' num2str(num_mod)]) == 1) && (contains(tmp(m).name, data_txt) == 1) ...
                    && (contains(tmp(m).name, vox_struc{voxels}) == 1)
                tmp_name = extractAfter(tmp(m).name, '_Fit_'); 
                true_RSquare = str2double(tmp_name(1:7));
        end
    end
    
    % One-tailed p-value
    p_perm = sum(RSquare_mat>=true_RSquare)/n_perm; % *2 if two-tailed
    fprintf('p-value: %f\t(true RSquare = %f)\n', p_perm, true_RSquare);
end

% Save results
t_day = datestr(now, 1);
t_time = strcat(datestr(now, 'HH'), '-', datestr(now, 'MM'));
if model_version ~= 3
    file_modav = [model_path '\ModelAveraging\ModelAveraging_' num2str(num_mod) '_' data_txt '_' vox_txt '_Fit_' num2str(RSquare_av) '_' t_day '_' t_time '.mat'];
    save(file_modav, 'mat_sort');
    file_meanpred = [model_path '\ModelAveraging\MeanPredictions_' data_txt '_' vox_txt '_' t_day '_' t_time '.txt'];
    dlmwrite(file_meanpred, [pred_mean test_sort]);
elseif model_version == 3
    file_modav = [model_path '\ModelAveraging\PermAveraging_' num2str(num_mod) '_' data_txt '_' vox_txt '_Fit_' num2str(perm_selected) '_p_' num2str(p_perm) '_' t_day '_' t_time '.mat'];
    save(file_modav, 'perm_sort_fit');
    file_meanpred = [model_path '\ModelAveraging\ModAv_perm_sort_fit_' data_txt '_' vox_txt '_' t_day '_' t_time '.txt'];
    dlmwrite(file_meanpred, perm_sort_fit);
end

%% Plot
if model_version ~= 3
    figure('Name', 'Model fit')
    hold on
    scatter(test_sort, pred_mean, 15, [0.4940, 0.1840, 0.5560], 'filled'); % Purple
    pol = polyfit(test_sort, pred_mean, 1);
    y1 = polyval(pol, test_sort);
    plot(test_sort, y1, 'LineWidth', 2); % Regression line
    text(0.42, 0.9, ['R^2 = ' num2str(RSquare_av)]);
    xlabel('Test outcome');
    ylabel('Mean predicted outcome');
    axis([0 1 0 1]);
    legend('Patient', 'Averaged model fit', 'Location', 'northwest');
    hold off
    saveas(gcf, strcat(model_path, '\ModelAveraging', '\Fig_ModelAveraging_', num2str(num_mod), '_', data_txt, '_', vox_txt, '_avRSquare_', num2str(RSquare_av), '_', t_day, '_', t_time, '.png'));    
end

