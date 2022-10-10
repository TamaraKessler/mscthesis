%% ******************* Data Analysis of WMH Using SVR *********************
%
% Lisa Röhrig, 12-02-2020, last update on 11-17-2021
%
% Applying Support Vector Regression for Prediction Analysis
% Can White Matter Hyperintensities predict neglect symptoms and severity  
% in stroke patients?
%
% ### MATRIX CONCATENATION [MERGING LESION + WMH MAPS] ###
%
% >> Models:
%       * Model 1 uses Lesion data 
%       * Model 2 uses Lesion + WMH data
%       * Model 3 uses Lesion + permuted WMH data
%       * Model 4 uses Lesion data + age
%
% >> WMH data formats:
%       * Type 1 uses volumetric data of the whole brain
%       * Type 2 uses volumetric data of the left hemisphere
%       * Type 3 uses volumetric data of the right hemisphere
%       * Type 4 uses rating scores (CHS/Fazekas) 
%
% >> Steps:
%       * Load data
%       * Permutation of data (if Model 3)
%       * Ignore non-lesioned voxels
%       * Principal component analysis
%       * Nested cross-validation
%       * Model fit calculation

clear; clc; close all;

% Path requirements for this script: within the main folder, there should be folders for
% structural data ('Data', with subfolders 'Hemisphere\WMH_left', 'Hemisphere\WMH_right',
% 'WholeBrain\Lesions' and 'WholeBrain\WMH') containing lesion maps in Nifti (.nii) format
% and for output ('Results'). The excel-file with WMH ratings, behavior etc. 
% should be within the 'Data'-folder. The 'Results'-folder should contain folders for
% each model version (e.g., 'Model_1', 'Model_2\Addition', 'Model_3\Concatenation'), 
% each of those subfolders should contain an additional 'Images'-folder.  
path = 'D:\Lisa\WMH\';
cd(path);

% Preparation for parfor (parallel computing)
delete(gcp('nocreate'))
try
    delete(myCluster.Jobs)
end

%% _____________________________ Load Data ________________________________
fprintf('MATRIX CONCATENATION\n');
fprintf('***SVR*** - Loading data...\n');

% Define mode
fprintf('Model versions: 1 - Lesion data, 2 - Lesion & WMH data, 3 - Lesion & permuted WMH data, 4 - Lesion data & age\n');
model_version = input('Select model version. Version # ');
while (model_version ~= 1) && (model_version ~= 2) && (model_version ~= 3) && (model_version ~= 4)
   model_version = input('Select model version. Number should be between 1 and 4! Version # ');
end
if (model_version == 1) || (model_version == 4)
    img_dir = fullfile(path, 'Results', ['Model_' num2str(model_version)], 'Images');
elseif (model_version == 2) || (model_version == 3)
    img_dir = fullfile(path, 'Results', ['Model_' num2str(model_version)], 'Concatenation', 'Images');
end

% Define data type
voxels = 0;
if (model_version ~= 1) && (model_version ~= 4) % Only if WMH data
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
    % Define which rating scale to use
    elseif data_format == 4
        fprintf('Rating scales: \t1 - CHS, 2 - Fazekas PV-WMH, 3 - Fazekas DS-WMH, 4 - PV-WMH & DS-WMH\n');
        rating = input('Select rating scale. Rating # ');
        while (rating ~= 1) && (rating ~= 2) && (rating ~= 3) && (rating ~= 4)
            rating = input('Select rating scale. Number should be between 1 and 4! # ');
        end
    end
elseif (model_version == 1) || (model_version == 4)
    data_format = 1;
end
switch data_format % For structure name later
    case 1
        if model_version == 1
            data_txt = 'VolWB';
        elseif (model_version == 2) || (model_version == 3)
            data_txt = 'Conc_VolWB';
        elseif model_version == 4
            data_txt = 'Conc_Age';
        end
    case 2
        data_txt = 'Conc_VolLHS';
    case 3
        data_txt = 'Conc_VolRHS';
    case 4
        if rating == 1
            data_txt = 'Conc_CHS'; 
        elseif rating == 2
            data_txt = 'Conc_PVWMH';
        elseif rating == 3
            data_txt = 'Conc_DSWMH';
        elseif rating == 4
            data_txt = 'Conc_Fazekas'; % PV-WMH and DS-WMH concatenated
        end
end
switch voxels 
    case 0
        vox_txt = [];
    case 1
        vox_txt = 'AllVox';
    case 2
        vox_txt = 'ShortVox';
end

% Random stream preparation
rs = input('Select number of random streams (usually between 1 and 10). rs = ');
while (abs(round(rs)-rs) == 0) == 0 % Only whole numbers
    rs = input('Select number of random streams (usually between 1 and 10). A whole number is needed! rs = ');
end

% Define number of folds during cross validation (last step of script)
k_oL = input('Choose k for a k-fold CV (OUTER LOOP; default is 10). k = ...');
if isempty(k_oL) == 1
    k_oL = 10; % Default
end
k_iL = input('Choose k for a k-fold CV (INNER LOOP; default is 5). k = ...');
if isempty(k_iL) == 1
    k_iL = 5; % Default
end
k = [k_oL, k_iL];

% Select number of permutations
if (model_version == 1) || (model_version == 2) || (model_version == 4)
    n_perm = 1;
elseif model_version == 3 % Permuted WMH data
    n_perm = input('\nPlease enter number of permutation-iterations (default = 5000). n = ...\n');
    if isempty(n_perm) == 1
        n_perm = 5000; % Default
    end
end

% .........................................................................
% Select excel-file containing data 
[num, txt, raw] = xlsread(fullfile(path, 'Data', 'WMH_Data'), 'WMH_Final');
if isequal(num(:, 1), sort(num(:, 1))) == 0
    fprintf('\nError! Data file has to be sorted regarding patient-number!!!\n\n');
    return
end

% Number of training examples (= number of patients)
m = size(num, 1); 

% Read normalized lesions/WMH 
% Table -> rows = patients, columns = predictors
for i = 1:m
    lesion_hdr = spm_vol(fullfile(path, 'Data', 'WholeBrain', 'Lesions', ...
            ['patient_' num2str(num(i, 1) + 1000) '_lesion.nii'])); % Patient-ID plus 1000 in name (1001, 1002...)
    lesion_img = spm_read_vols(lesion_hdr);
    if (model_version == 2) || (model_version == 3) % WMH data
        if data_format == 1 % Volumetric data of whole brain
            wmh_hdr = spm_vol(fullfile(path, 'Data', 'WholeBrain', 'WMH', ...
                ['patient_' num2str(num(i, 1) + 1000) '_WMH.nii']));
            wmh_img = spm_read_vols(wmh_hdr);
        elseif data_format == 2 % Volumetric data of left HS
            wmh_hdr = spm_vol(fullfile(path, 'Data', 'Hemisphere', 'WMH_left', ...
                ['patient_' num2str(num(i, 1) + 1000) '_WMH_leftHS.nii']));
            wmh_img = spm_read_vols(wmh_hdr);
        elseif data_format == 3 % Volumetric data of right HS
            wmh_hdr = spm_vol(fullfile(path, 'Data', 'Hemisphere', 'WMH_right', ...
                ['patient_' num2str(num(i, 1) + 1000) '_WMH_rightHS.nii']));
            wmh_img = spm_read_vols(wmh_hdr);
        end
    end
    
    % Vectorize 3D-matrix of lesion/WMH maps for each patient
    lesion_v = uint8(lesion_img(:)');
    lesions_all(i, :) = lesion_v;
    if ((model_version == 2) || (model_version == 3)) ... % WMH data
            && (data_format ~= 4)
        wmh_v = uint8(wmh_img(:)');
        wmh_org(i, :) = wmh_v;
    end
    
    clearvars lesion_img wmh_img lesion_hdr lesion_v wmh_hdr wmh_v
end

% Define output variables
outcome_org = num(:, 20); % Mean CoC 
outcome_pos = outcome_org; 
outcome_pos(outcome_pos < 0) = 0; % Only positive CoC (negative ones -> 0)
outcome = realsqrt(outcome_pos);

% Age
if model_version == 4
    age = num(:, 2);
    % Feature scaling
    age_min = min(age, [], 'all');
    age_max = max(age, [], 'all');
    age_org = (age - age_min)/(age_max - age_min);
    fprintf('Age range after feature scaling: Min = %d, Max = %d\n', ...
        min(age_org, [], 'all'), max(age_org, [], 'all'))
end

% Define scale values
if exist('rating') == 1
    if rating == 1
        wmh_rat = num(:, 12); % CHS
    elseif rating == 2
        wmh_rat = num(:, 13); % PV-WMH
    elseif rating == 3
        wmh_rat = num(:, 14); % DS-WMH
    elseif rating == 4
        wmh_rat = [num(:, 13) num(:, 14)]; % PV-WMH & DS-WMH concatenated
    end
    
    % Feature scaling (values between 0 and 9)
    sc_min = min(wmh_rat, [], 'all');
    sc_max = max(wmh_rat, [], 'all');
    wmh_org = (wmh_rat - sc_min)/(sc_max - sc_min);
    fprintf('Rating scores (WMH extent) range after feature scaling: Min = %d, Max = %d\n', ...
        min(wmh_org, [], 'all'), max(wmh_org, [], 'all'));
end

lesions_org = lesions_all;

%%
% For all iterations
tic
myPool = parpool(4); % If an error "out of memory" occurs, the number of pools must be reduced.
parfor r = 1:rs
% Preparation for parfor-loop
perm_list = []; perm_pred = []; pca_struc = [];
fprintf('\n********** PARFOR-ITERATION # %d **********\n', r);

for p = 1:n_perm
    %% ________________________ Permutation of Data _______________________
    % Begin of new permutation iteration
    lesions_all = lesions_org;
    if (model_version ~= 1) && (model_version ~= 4)
        wmh_all = wmh_org;
    elseif model_version == 4
        age_all = age_org;
    end
    if model_version == 3
        wmh_short = []; wmh_s = []; wmh_s_log = []; wmh_ind = [];

        % Permuted WMH data
        defaultStream = RandStream('mlfg6331_64');
        defaultStream.Substream = p;
        [wmh_perm, perm_idx] = permutation(wmh_org, defaultStream); %disp(perm_idx);
        fprintf('Permutation of WMH data (iteration # %d)..... done\n', p);
        wmh_all = wmh_perm;
        defaultStream = [];
    end
    
    %% _____________________ Ignore Non-Lesioned Voxels _______________________
    fprintf('***SVR*** - Keep only lesioned voxels...\n');

    % Sum of patients that have lesion in specific voxel
    lesions_s = sum(lesions_all, 1);
    lesions_s_log = (lesions_s > 4); % Voxel is lesioned in at least 5 patients!!
    if ((model_version == 2) || (model_version == 3)) && (data_format ~= 4)
            wmh_s = sum(wmh_all, 1);
        if voxels == 1 % Keep voxels that are damaged in at least 1 patient
            wmh_s_log = (wmh_s > 0);
        elseif voxels == 2 % Keep voxels that are damaged in at least 5 patients
            wmh_s_log = (wmh_s > 4);
        end
    end

    % Only keep columns/voxels that are of interest
    % To keep variance small, select voxels for lesion and WMH map separately
    lesions_short = lesions_all(:, lesions_s_log == 1);
    if ((model_version == 2) || (model_version == 3)) && (data_format ~= 4) % WMH
        wmh_short = wmh_all(:, wmh_s_log == 1);
    end
    
    % Put all maps into one matrix
    if ((model_version == 2) || (model_version == 3)) && (data_format ~= 4) % WMH data
        maps_short = [lesions_short wmh_short];
    end
    
    lesions_s = []; lesions_s_log = []; wmh_s = []; wmh_s_log = []; wmh_short = [];
    if model_version ~= 3
        lesions_all = []; 
        if data_format ~= 4
            wmh_all = [];
        end
    end
    
    %% ___________________ Principal Component Analysis __________________
    fprintf('***** Principal Component Analysis *****\n');
    pca_log = 1; % Change to 0 if PCA is not desired

    % Feature reduction to save disk space/computation effort 
    % (Preparation for permutation loops of nested CVs)
    if pca_log == 1
        if (model_version == 1) || (data_format == 4) || (model_version == 4)
            [coeff, score, latent, ~, explained] = pca(double(lesions_short));
            n_pca = size(score, 2); % Total number of principal components 
            %plot(sqrt(sum(score.^2, 1))); % Visual inspection
            old_lesions_short = lesions_short;
        elseif ((model_version == 2) || (model_version == 3)) && (data_format ~= 4)
            [coeff, score, latent, ~, explained] = pca(double(maps_short));
            n_pca = size(score, 2); % Total number of principal components
            %plot(sqrt(sum(score.^2, 1))); % Visual inspection
            old_maps_short = maps_short;
        end
        % Number of principal components that explain 98% of variance
        var_thr = 98; % Variance threshold
        cum_var = cumsum(explained); % Cumulative sum of variance explained by components
        for h = 1:length(cum_var)
            if cum_var(h) > var_thr
                n_comp = h; % Number of components
                fprintf('\nMaintain %d components (total = %d components) with %f cumulative percent explained variance\n\n', n_comp, n_pca, cum_var(h));
                break
            end
        end
        % Only keep components that explain at least 98% of variance
        if (model_version == 1) || (data_format == 4) || (model_version == 4)
            lesions_short = score(:, 1:n_comp); % Features in principal component space
        elseif ((model_version == 2) || (model_version == 3)) && (data_format ~= 4)
            maps_short = score(:, 1:n_comp); % Features in principal component space
        end

        % Make plot
        if (model_version == 1) || (data_format == 4) || (model_version == 4)
            title = 'Lesion data';
        elseif (model_version == 2) && (data_format ~= 4)
            title = 'Lesion + WMH data';
        end
        t_day = datestr(now, 1);
        t_time = strcat(datestr(now, 'HH'), '-', datestr(now, 'MM'));
        if ((model_version == 1) || (model_version == 2)) && (data_format ~= 4) 
            figure('Name', title)
            l = length(cum_var);
            line(1:l, cum_var, 'marker', 's', 'color', 'k', 'markerfacecolor', 'white', 'MarkerSize', 5); 
            yline(var_thr, '--k'); % Add threshold line
            xlabel('Number of principal components');
            ylabel('Cumulative variance explained');
            legend('Data', 'Threshold', 'Location', 'southeast');
            saveas(gcf, strcat(img_dir, '\Fig_PCA_', title, '_rs_', num2str(r), '_', data_txt, '_', t_day, '_', t_time, '.png')); 
        end

        % Save structure
        if ((model_version == 1) || (model_version == 2)) && (data_format ~= 4)
            pca_struc.model = model_version;
            if model_version == 1
                pca_struc.origDim = [num2str(size(old_lesions_short, 1)) '*' ...
                    num2str(size(old_lesions_short, 2))];
                pca_struc.redDim = [num2str(size(lesions_short, 1)) '*' ...
                    num2str(size(lesions_short, 2))];
            elseif (model_version == 2) 
                pca_struc.origDim = [num2str(size(old_maps_short, 1)) '*' ...
                    num2str(size(old_maps_short, 2))];   
                pca_struc.redDim = [num2str(size(maps_short, 1)) '*' ...
                    num2str(size(maps_short, 2))];
            end
            pca_struc.totalPCs = n_pca;
            pca_struc.redPCs = n_comp;
            pca_struc.explVar = cum_var(h);
            if model_version == 1
                file_pca = [path 'Results\Model_' num2str(model_version) '\PCA_struc_rs_' num2str(r) '_' data_txt '_' vox_txt '_' t_day '_' t_time '.mat'];
            elseif model_version == 2
                file_pca = [path 'Results\Model_' num2str(model_version) '\Concatenation\PCA_struc_rs_' num2str(r) '_' data_txt '_' vox_txt '_' t_day '_' t_time '.mat'];
            end
            file_pca_mat = matfile(file_pca, 'writable', true);
            file_pca_mat.pca_struc = pca_struc;
        end

        % Feature Scaling (needed, because features have a large range
        % after PCA; desired range: 0 to 1) >> Min-Max-Normalization
        if (model_version == 1) || (data_format == 4) || (model_version == 4) % Lesion
            les_min = min(lesions_short, [], 'all');
            les_max = max(lesions_short, [], 'all');
            les_scaled = (lesions_short - les_min)/(les_max - les_min);
            lesions_short = les_scaled;
            fprintf('Feature range after scaling: Min = %d, Max = %d\n', ...
                min(lesions_short, [], 'all'), max(lesions_short, [], 'all'));
        elseif ((model_version == 2) || (model_version == 3)) && (data_format ~= 4) % WMH
            maps_min = min(maps_short, [], 'all');
            maps_max = max(maps_short, [], 'all');
            maps_scaled = (maps_short - maps_min)/(maps_max - maps_min);
            maps_short = maps_scaled;
            fprintf('Feature range after scaling: Min = %d, Max = %d\n', ...
                min(maps_short, [], 'all'), max(maps_short, [], 'all'));
        end           
    end
    
    % Add WMH rating scores or age values to lesion-matrix
    if data_format == 4
       maps_short = [double(lesions_short) wmh_all];
    elseif model_version == 4
        maps_short = [double(lesions_short) age_all];
    end
    
    %% _____________________ Nested Cross Validation __________________________
    fprintf('***** Nested Cross-Validation *****\n');

    log2C = -5:1:15; % Ascending to pick the biggest value before reaching plateau
    log2G = -15:1:5; 
    
    % Prepare order of subjects
    % Define seed for random number selection
    defaultStream = RandStream('mlfg6331_64');
    if model_version == 1
        subjects_ind = (1:size(lesions_short, 1))';
        defaultStream.Substream = r;
        subjects_rand = permutation(subjects_ind, defaultStream);
        lesions_rand = lesions_short(subjects_rand, :);
    elseif (model_version == 2) || (model_version == 3) || (model_version == 4)
        subjects_ind = (1:size(maps_short, 1))';
        defaultStream.Substream = r;
        subjects_rand = permutation(subjects_ind, defaultStream);
        maps_rand = maps_short(subjects_rand, :);
    end
    outcome_rand = outcome(subjects_rand, :);
    defaultStream = [];

    % Nested CV - k-fold outer loop CV and (k/2)-fold inner loop CV
    % Minimizes MAE (mean absolute error) in inner loop and outputs RSquare (= Coefficient 
    % of Determination) in outer loop
    if model_version == 1 % Lesion data
       model_struc = nested_CV(lesions_rand, outcome_rand, subjects_rand, k, log2C, log2G, model_version);
    elseif (model_version == 2) || (model_version == 3) || (model_version == 4) % WMH
       model_struc = nested_CV(maps_rand, outcome_rand, subjects_rand, k, log2C, log2G, model_version);        
    end
    
    %% _______________________________ Model Fit ______________________________
    % Save all predictions and test-outcome in lists (all patients)
    pred_all = [];
    test_all = [];
    subj_all = [];
    for s = 1:size(model_struc, 2)
        pred_all = [pred_all; model_struc(s).predicted_label_oL];
        test_all = [test_all; model_struc(s).outcome_test_oL];
        subj_all = [subj_all; model_struc(s).subjects_test_oL];
    end

    % Calculate RSquare average
    squRes = (test_all - pred_all).^2;
    squTotal = (test_all - mean(test_all)).^2;
    RSquare_all = 1 - (sum(squRes)/sum(squTotal));
    if (model_version == 3) && (p == 1)
        diary([path 'Results\Model_' num2str(model_version) '\Concatenation\Permutation_output_' num2str(n_perm) '_' data_txt '_' vox_txt '_rs_' num2str(r) '_' t_day '_' t_time '.txt']);
    elseif (model_version == 3) && (p ~= 1)
        diary ON
    end
    if model_version == 3
        fprintf('\nRandomStream #%d; Permutation iteration #%d; RSquare = %g\n\n', r, p, RSquare_all);
        diary OFF
    else
        fprintf('\nRandomStream #%d; RSquare = %g\n\n', r, RSquare_all);
    end

    % Save plot showing coefficient of determination (predicted vs. real outcome)
    t_day = datestr(now, 1);
    t_time = strcat(datestr(now, 'HH'), '-', datestr(now, 'MM'));
    if model_version ~= 3 % Only if just one model fitted
        figure('Name', 'Model Fit')
        hold on
        scatter(test_all, pred_all, 15, [0.4940, 0.1840, 0.5560], 'filled'); % Purple
        pol = polyfit(test_all, pred_all, 1);
        y1 = polyval(pol, test_all);
        plot(test_all, y1, 'LineWidth', 2); % Regression line (unsure, whether correct one is displayed)
        text(0.35, 0.9, ['R^2 = ' num2str(RSquare_all)]);
        xlabel('Real outcome');
        ylabel('Predicted outcome');
        axis([0 1 0 1]);
        legend('Patient', 'Model Fit', 'Location', 'northwest');
        hold off
        saveas(gcf, strcat(img_dir, '\Fig_ModelFit_', data_txt, '_', vox_txt, '_rs_', num2str(r), '_RSquare_', num2str(RSquare_all), '_', t_day, '_', t_time, '.png'));    

        % Save structure
        if (model_version == 1) || (model_version == 4)
            file_nestedCV = [path 'Results\Model_' num2str(model_version) '\model_struc_nestedCV_' data_txt '_' vox_txt '_rs_' num2str(r) '_avRSquare_' num2str(RSquare_all) '_' t_day '_' t_time '.mat'];
        elseif model_version == 2
            file_nestedCV = [path 'Results\Model_' num2str(model_version) '\Concatenation\model_struc_nestedCV_' data_txt '_' vox_txt '_rs_' num2str(r) '_avRSquare_' num2str(RSquare_all) '_' t_day '_' t_time '.mat'];
        end
        model_struc_mat = matfile(file_nestedCV, 'writable', true);
        model_struc_mat.model_struc = model_struc;
    end
    
    % Permutation testing
    if model_version == 3
        % Save list of model fits
        perm_list(p, 1) = p;
        perm_list(p, 2) = n_comp;
        perm_list(p, 3) = RSquare_all;
        % Save predictions
        [subj_sort, sort_ind] = sort(subj_all);
        pred_sort = pred_all(sort_ind, :);
        test_sort = test_all(sort_ind, :);        
        perm_pred(:, p) = pred_sort;
        % Close all figures
        close all
    end
end

%%
% Permutation results
if model_version == 3
    % Sort list 
    [perm_sort_fit, perm_ind] = sort(perm_list(:, 3), 1, 'descend');
    perm_sort_other = perm_list(perm_ind, 1:2);
    perm_sort = [perm_sort_other perm_sort_fit];
    % Selecting 95%-best model
    num_p = length(perm_sort) - (length(perm_sort) * 0.95);
    perm_selected = perm_sort(ceil(num_p), :);
    % Save list 
    file_perm = [path 'Results\Model_' num2str(model_version) '\Concatenation\Permutation_results_' num2str(p) '_' data_txt '_' vox_txt '_Fit_' num2str(perm_selected(1, 3)) '_rs_' num2str(r) '_' t_day '_' t_time '.mat'];
    perm_sort_mat = matfile(file_perm, 'writable', true);
    perm_sort_mat.perm_sort = perm_sort;    
    % Save prediction matrix
    file_pred = [path 'Results\Model_' num2str(model_version) '\Concatenation\Permutation_predictions_' num2str(n_perm) '_' data_txt '_' vox_txt '_Fit_' num2str(perm_selected(1, 3)) '_rs_' num2str(r) '_' t_day '_' t_time '.mat'];
    file_pred_mat = matfile(file_pred, 'writable', true);
    file_pred_mat.perm_pred = perm_pred;
end
end
fprintf('\nEnd of script.\n\n');
toc
