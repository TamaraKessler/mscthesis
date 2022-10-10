function [model_struc] = nested_CV(features_rand, outcome_rand, subjects_rand, k, log2C, log2G, model_version)
% *********************** Nested Cross Validation *************************
%
% Lisa Röhrig, 12-16-2020, last update on 06-08-2021
%
% Nested Cross-Validation (k-fold)
% Inner loop:   model fitting procedure (grid search)
% Outer loop:   estimates performance of model created by inner loop,
%               splitting data into training and test-sets
% Model selection:  * take C and gamma-values used in the best
%                   performing model (inner loop), regardless of calculated 
%                   performance (outer loop)!
%                   * if estimated performance of inner loop and outer loop
%                   differs a lot (best model), then model might be
%                   overfitted (inner loop)
% Input parameters:     features - data set X
%                       outcome - behavioural data Y
%                       k - number of folds used during cross validation (k-fold)
%                       [k(1) - k-value outer loop, k(2) - k-value inner loop]                      
%                       log2C - range of values to search for best C
%                       log2G - range of values to search for best gamma
%                       model_version - defines type of model (used data
%                           type, permutation testing, ...)
% Output parameters:    model_struc - structure containing relevant
%                           variables regarding inner and outer loop models
%
% Function uses libsvm package by Chang & Lin (2011), available at
% http://www.csie.ntu.edu.tw/~cjlin/libsvm).

% Define group sizes
split_numb = size(features_rand, 1) / k(1);
group_size_mean = floor(split_numb);
rest_numb = size(features_rand, 1) - (k(1) * group_size_mean);
group_sizes_v = zeros(k(1), 1);
group_sizes_v(1:(k(1)-rest_numb), 1) = group_size_mean;
group_sizes_v((k(1)-rest_numb+1):end, 1) = group_size_mean + 1;
if sum(group_sizes_v) ~= size(features_rand, 1) % Check
    fprintf('\nERROR. Wrong group sizes. \n');
    return
end

% Outer loop - estimation of model performance
for ii = 1:k(1)
    fprintf('\n *** Nested cross-validation... outer loop iteration # %d\n\n', ii);
    
    % Split data into k-folds
    group_size_test = group_sizes_v(ii, 1);
    group_size_train = sum(group_sizes_v)-group_sizes_v(ii, 1);
    features_test = features_rand(1 : group_size_test, :); % 1 part serves as test set
    outcome_test = outcome_rand(1 : group_size_test, :);
    subjects_test = subjects_rand(1 : group_size_test, :);
    features_train = features_rand(group_size_test+1 : end, :); % k-1 parts merged to one training set
    outcome_train = outcome_rand(group_size_test+1 : end, :);
    subjects_train = subjects_rand(group_size_test+1 : end, :);
    
    % Move current test-set to end of matrix (preparation for next
    % iteration)
    features_rand = [features_rand(group_size_test+1 : end, :); ...
        features_rand(1:group_size_test, :)];
    outcome_rand = [outcome_rand(group_size_test+1 : end, :); ...
        outcome_rand(1:group_size_test, :)];
    subjects_rand = [subjects_rand(group_size_test+1 : end, :); ...
        subjects_rand(1:group_size_test, :)];
    
    % Inner loop - model selection / course grid search
    % "Train-set" of outer loop is split into train and test sets for inner
    % loop during grid search
    if model_version == 3 % Silent mode
        [perform_list, params_list, cv_struc, bestC, bestG, bestMAE, bestCorr] = GridSearch_coeffDeterm(outcome_train, features_train, subjects_train, log2C, log2G, k(2), 'off'); 
    else
        [perform_list, params_list, cv_struc, bestC, bestG, bestMAE, bestCorr] = GridSearch_coeffDeterm(outcome_train, features_train, subjects_train, log2C, log2G, k(2), 'on'); 
    end
        
    % Performance calculation of model selected during inner loop
    options = ['-s 3 -t 2 -p 0.1 -c ' num2str(bestC) ' -g ' num2str(bestG) ' -q']; % Without cv; '-q' to disable screen output    
    model = svmtrain(outcome_train, sparse(double(features_train)), options); % Need structure for prediction
    [results, predicted_label, accuracy, ~] = evalc('svmpredict(outcome_test, sparse(double(features_test)), model)');

    i_start = strfind(results, 'Mean squared error = ') + length('Mean squared error = ');
    i_end = strfind(results, ' (regression)') - 1; % Next line begins with 'Squared correlation coefficient'
    MSE_oL = str2double(results(i_start : i_end(1)));
    MAE_oL = sqrt(MSE_oL); 
        
    % Save parameters in structure
    model_struc(ii).iteration_oL = ii;
    model_struc(ii).subjects_test_oL = subjects_test;
    model_struc(ii).bestC_iL = bestC;
    model_struc(ii).bestG_iL = bestG;
    model_struc(ii).bestMAE_iL = bestMAE;
    model_struc(ii).model_MAE_oL = MAE_oL;
    model_struc(ii).outcome_test_oL = outcome_test;
    model_struc(ii).predicted_label_oL = predicted_label;
    if model_version ~= 3 % Only if not permutation
        model_struc(ii).group_sizes_train_oL = group_size_train;
        model_struc(ii).group_size_test_oL = group_size_test;
        model_struc(ii).perform_list_iL = perform_list;
        model_struc(ii).params_list_iL = params_list;
        model_struc(ii).cv_iL = cv_struc;
        model_struc(ii).bestCorr_iL = bestCorr;
        model_struc(ii).model_MSE_oL = MSE_oL; 
        model_struc(ii).model_Corr_oL = accuracy(3);
        model_struc(ii).accuracy_oL = accuracy;
    
        fprintf('\nModel # %d:\ninner loop -> best C = %g, best gamma = %g, best MAE = %g\nouter loop -> Corr = %f\n', ...
            ii, bestC, bestG, bestMAE, accuracy(3));
    end
end
end
