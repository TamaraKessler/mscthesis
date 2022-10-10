function [perform_list, params_list, cv_struc, bestC, bestG, bestMAE, bestCorr] = GridSearch_coeffDeterm(outcome, features, subjects, log2C, log2G, k, mode)
% Function GridSearch_coeffDeterm optimizes hyperparameters
% with k-fold CV on different values and combinations of cost parameter (C) 
% and gamma in epsilon-SVR machine learning algorithms with RBF kernel.
% 
% Input argument: outcome - column vector of outcome variable
%                 features - scaled feature matrix (observations x
%                   features)
%                 log2C - vector of C-values used in parameter
%                   optimization, format: log2(C)
%                 log2G - vector of gamma-values used in paramter
%                   optimization, format: log2(gamma)
%                 k - optional, k-fold cross validation (default is k = 5)
%                 mode - 'on' with text messages in command window, 'off'
%                   silent
%
% Output argument: perform - matrix for plotting of dimensions C x gamma
%                  params - parameter matrix with columns C, gamma, R^2 
%                  bestC - best parameter C
%                  bestG - best parameter gamma
%                  bestMAE - lowest mean absolute error (best model fit)
%
% Original version written by Vanessa Kasties on Sept 28 2020. 
% Modified by Lisa Röhrig (Dec 2020, last updated on 25-01-2021):  
% Function minimizes mean absolute error (MAE).
% Function uses libsvm package by Chang & Lin (2011), available at
% http://www.csie.ntu.edu.tw/~cjlin/libsvm).

num_C = length(unique(log2C)); % Number of unique Cs
num_G = length(unique(log2G)); % Number of unique gammas

% Degree of cross validation (k-fold)
if exist('k', 'var') == 0
    k = 5; % Default
end

% Initialize best prediciton accuracy
bestMAE = 1000; % High initial MAE, function minimizes error
bestC = 0;
bestG = 0;
bestCorr = 0;
i = 1;

% Initialize matrix for plotting Rsquare for each C-gamma pair; C rows, gamma columns
perform_list.MSE = zeros(num_C, num_G); 
perform_list.Corr = zeros(num_C, num_G);

% Initialize matrix to save parameters and Rsquare
params_list = zeros(num_C*num_G, 3); 

% Define group sizes
split_numb = size(features, 1) / k;
group_size_mean = floor(split_numb);
rest_numb = size(features, 1) - (k * group_size_mean);
group_sizes_v = zeros(k, 1);
group_sizes_v(1:(k-rest_numb), 1) = group_size_mean;
group_sizes_v((k-rest_numb+1):end, 1) = group_size_mean + 1;
if sum(group_sizes_v) ~= size(features, 1) % Check
    fprintf('\nERROR. Wrong group sizes. \n');
    return
end

% Inner loop of nested CV (each parameter mix is tested on all cv-folds)
for pc = 1:length(log2C)
    C_item = log2C(pc);
    for pg = 1:length(log2G)
        G_item = log2G(pg);

        % Initialize ML parameters
        cmd = ['-s 3 -t 2 -g ', num2str(2^G_item), ' -c ', num2str(2^C_item), ' -p 0.1 -q']; % Adjust p for epsilon-value (margin size)

        % For each fold
        for ii = 1:k
            % Select fold
            group_size_test = group_sizes_v(ii, 1);
            group_size_train = sum(group_sizes_v)-group_sizes_v(ii, 1);
            features_test = features(1 : group_size_test, :); % 1 part serves as test set
            outcome_test = outcome(1 : group_size_test, :);
            subjects_test = subjects(1 : group_size_test, :);
            features_train = features(group_size_test+1 : end, :); % k-1 parts merged to one training set
            outcome_train = outcome(group_size_test+1 : end, :);
            subjects_train = subjects(group_size_test+1 : end, :);

            % Move current test-set to end of matrix (preparation for next iteration)
            features = [features(group_size_test+1 : end, :); ...
                features(1:group_size_test, :)];
            outcome = [outcome(group_size_test+1 : end, :); ...
                outcome(1:group_size_test, :)];
            subjects = [subjects(group_size_test+1 : end, :); ...
                subjects(1:group_size_test, :)];

            % Train model
            model = svmtrain(outcome_train, sparse(double(features_train)), cmd); 
            
            % Predict outcome (Saves output to command window in results)
            [results, predicted_label, accuracy, ~] = evalc('svmpredict(outcome_test, sparse(double(features_test)), model)');        
            i_start1 = strfind(results, 'Squared correlation coefficient = ') + length('Squared correlation coefficient = ');
            i_end = strfind(results, ' (regression)') - 1;
            RSquare_Corr = str2double(results(i_start1 : i_end(2)));
            i_start2 = strfind(results, 'Mean squared error = ') + length('Mean squared error = ');
            MSE = str2double(results(i_start2 : i_end(1))); % (outcome-predicted).^2
            MAE = sqrt(MSE); % (outcome-predicted)
            
            % Save results of current fold
            cv_fold(ii).MSE = MSE;
            cv_fold(ii).MAE = MAE;
            cv_fold(ii).Corr = RSquare_Corr;
            if strcmp(mode, 'on') == 1
                cv_fold(ii).train_subj = subjects_train;
                cv_fold(ii).test_out = sparse(outcome_test);
                cv_fold(ii).test_subj = subjects_test;
                cv_fold(ii).predicted = sparse(predicted_label);
                cv_fold(ii).accuracy = accuracy;
            end
        end

        % Calculate mean error/correlation of cv
        MAE_cv = mean([cv_fold.MAE]);
        MSE_cv = mean([cv_fold.MSE]);
        Corr_cv = mean([cv_fold.Corr]);
        
        % Check if current hyperparameters revealed best accuracy
        if (MAE_cv <= bestMAE) % Minimize mean absolute error (MAE)
            bestMAE = MAE_cv;
            %bestMSE = MSE_cv;
            bestC = 2^C_item;
            bestG = 2^G_item;
            bestCorr = Corr_cv;
        end
        
        % Print mean cv results for current hyperparameter combination
        if strcmp(mode, 'on') == 1
            fprintf('%g %g %g %g (best C = %g, gamma = %g, MAE = %g, R²-Corr = %g)\n', C_item, G_item, MAE, RSquare_Corr, bestC, bestG, bestMAE, bestCorr);
        end
        perform_list.MAE(pc, pg) = MAE_cv; % loop1: C1, g1, loop2: C1, g2, ...
        perform_list.Corr(pc, pg) = Corr_cv;
        params_list(pg + (pc-1)*num_G, 1) = i;
        params_list(pg + (pc-1)*num_G, 2) = C_item;
        params_list(pg + (pc-1)*num_G, 3) = G_item;
        params_list(pg + (pc-1)*num_G, 4) = Corr_cv;
        params_list(pg + (pc-1)*num_G, 5) = MSE_cv;
        params_list(pg + (pc-1)*num_G, 6) = MAE_cv;
        cv_struc(i).cv_folds = cv_fold;
        cv_struc(i).C_item = C_item;
        cv_struc(i).G_item = G_item;

        i = i+1;
    end
end

% Check for good balance between small error (MAE) and high correlation
% (RSquare_Corr)
[params_sort, sort_idx] = sort(params_list(:, 6)); % Sort MAE
params2 = params_list(:, 1:5);
params2 = params2(sort_idx, :);
params_list = [params2 params_sort];

% If multiple models result in smallest error, select best model with 
% smallest error and highest correlation
n_sMAE = sum(params_list(:, 6) == params_list(1, 6)); % Amount of smallest error
[max_Corr, max_idx] = max(params_list(1:n_sMAE, 4)); % Looking for greatest correlation

% Best model (best error-correlation balance)
bestC = 2^params_list(max_idx, 2);
bestG = 2^params_list(max_idx, 3);
bestCorr = max_Corr;
bestMAE = params_list(max_idx, 6);
fprintf('\n***** Model with best error-correlation-balance selected:');
fprintf('\nC = %g (%g), gamma = %g (%g), MAE = %g (MSE = %g), R²-Corr = %g *****\n', ...
    bestC, params_list(max_idx, 2), bestG, params_list(max_idx, 3), bestMAE, params_list(max_idx, 5), bestCorr);
end