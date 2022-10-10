
%% Set Up

clear all; clc;
addpath('D:\MATLAB\libsvm-3.24\libsvm-3.24\matlab');
addpath('C:\Users\Sperber\Desktop\Sex Diff Scripts\OGs');
datadir = 'D:\Tamara\Diskonnektionen\LibSVM\lesionmaps';
cd(datadir);

instmat = readmatrix(fullfile(datadir, 'instance_matrix.csv'));
sparse_instmat = sparse(instmat);

[N, M] = size(sparse_instmat);

labels = readmatrix(fullfile(datadir, 'labels_4class.csv'));

%% Pre-Allocate some memory

defaultStream = RandStream('mlfg6331_64');
% log2C = -5:1:14; % Ascending to pick the biggest value before reaching plateau
% Nu = 0.01:0.05:1;

allsubs = 1:N;
allsubs = allsubs';

oL = 10;
iL = 5;

log2C = -5:1:15;
Nu = 0.01:0.05:0.51;

num_C = length(unique(log2C)); % Number of unique Cs
num_Nu = length(unique(Nu)); % Number of unique nus

cv_fold = struct;
cv_struc = struct;

params_list = zeros(num_C*num_Nu, 3); 

model_struc = struct;

%% Define group sizes

idx_data = (1:size(sparse_instmat, 1))';
idx_perm = randperm(defaultStream, length(idx_data))';
features_rand = sparse_instmat(idx_perm, :);
outcome_rand = labels(idx_perm, :);
subjects_rand = allsubs(idx_perm, :);
    
split_numb = size(features_rand, 1) / oL;
group_size_mean = floor(split_numb);
rest_numb = size(features_rand, 1) - (oL(1) * group_size_mean);
group_sizes_v = zeros(oL(1), 1);
group_sizes_v(1:(oL(1)-rest_numb), 1) = group_size_mean;
group_sizes_v((oL(1)-rest_numb+1):end, 1) = group_size_mean + 1;
if sum(group_sizes_v) ~= size(features_rand, 1) % Check
    fprintf('\nERROR. Wrong group sizes. \n');
    return
end

%%

for i_oL = 1:oL
    
    % Split data into training & testing set
    group_size_test = group_sizes_v(i_oL, 1);
    group_size_train = sum(group_sizes_v)-group_sizes_v(i_oL, 1);
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
    
    %% inner loop
    
    best_acc = 30;
    i = 1;
    
    tmp = struct;
    tmp.outcome = outcome_train;
    tmp.features = features_train;
    tmp.subjects = subjects_train;
    
    iL_split_numb = size(tmp.features, 1) / iL;
    iL_group_size_mean = floor(iL_split_numb);
    iL_rest_numb = size(tmp.features, 1) - (iL * iL_group_size_mean);
    iL_group_sizes_v = zeros(iL, 1);
    iL_group_sizes_v(1:(iL-iL_rest_numb), 1) = iL_group_size_mean;
    iL_group_sizes_v((iL-iL_rest_numb+1):end, 1) = iL_group_size_mean + 1;
    if sum(iL_group_sizes_v) ~= size(tmp.features, 1) % Check
        fprintf('\nERROR. Wrong group sizes. \n');
        return
    end
    
    % Inner loop of nested CV (each parameter mix is tested on all cv-folds)
    for pc = 1:length(log2C)
        C_item = log2C(pc);
        for pnu = 1:length(Nu)
            Nu_item = Nu(pnu);
            
            % Initialize ML parameters
            cmd = ['-s 1 -t 2 -n ', num2str(Nu_item), ' -c ', num2str(2^C_item), ' -q']; % Adjust p for epsilon-value (margin size)
            
            % For each fold
            for i_iL = 1:iL
                % Select fold
                iL_group_size_test = iL_group_sizes_v(i_iL, 1);
                iL_group_size_train = sum(iL_group_sizes_v)-iL_group_sizes_v(i_iL, 1);
                iL_features_test = tmp.features(1 : iL_group_size_test, :); % 1 part serves as test set
                iL_outcome_test = tmp.outcome(1 : iL_group_size_test, :);
                iL_subjects_test = tmp.subjects(1 : iL_group_size_test, :);
                iL_features_train = tmp.features(iL_group_size_test+1 : end, :); % k-1 parts merged to one training set
                iL_outcome_train = tmp.outcome(iL_group_size_test+1 : end, :);
                iL_subjects_train = tmp.subjects(iL_group_size_test+1 : end, :);
                
                % Move current test-set to end of matrix (preparation for next iteration)
                tmp.features = [tmp.features(iL_group_size_test+1 : end, :); ...
                    tmp.features(1:iL_group_size_test, :)];
                tmp.outcome = [tmp.outcome(iL_group_size_test+1 : end, :); ...
                    tmp.outcome(1:iL_group_size_test, :)];
                tmp.subjects = [tmp.subjects(iL_group_size_test+1 : end, :); ...
                    tmp.subjects(1:iL_group_size_test, :)];
                
                % Train model
                model = svmtrain(iL_outcome_train, sparse(double(iL_features_train)), cmd);
                
                % Predict outcome (Saves output to command window in results)
                [results, predicted_label, accuracy, ~] = evalc('svmpredict(iL_outcome_test, sparse(double(iL_features_test)), model)');
                num_correct = extractBetween(results, '% (', ') (classification)');
                
                % Save results of current fold
                cv_fold(i_iL).accuracy = accuracy(1);
                cv_fold(i_iL).number_correct = num_correct{:};
                cv_fold(i_iL).train_subj = subjects_train;
                cv_fold(i_iL).test_out = sparse(iL_outcome_test);
                cv_fold(i_iL).test_subj = iL_subjects_test;
                cv_fold(i_iL).predicted = sparse(predicted_label);
            end
            
            if accuracy(1) >= best_acc
                best_acc = accuracy(1);
                bestC = 2^C_item;
                bestNu = Nu_item;
            end
          
            % Print mean cv results for current hyperparameter combination
            %fprintf('%g %g %g %g (best C = %g, Nu = %g)\n', C_item, Nu_item, bestC, bestNu);
            
            params_list(pnu + (pc-1)*num_Nu, 1) = i;
            params_list(pnu + (pc-1)*num_Nu, 2) = C_item;
            params_list(pnu + (pc-1)*num_Nu, 3) = Nu_item;
            params_list(pnu + (pc-1)*num_Nu, 4) = accuracy(1);
            cv_struc(i).cv_folds = cv_fold;
            cv_struc(i).C_item = C_item;
            cv_struc(i).Nu_item = Nu_item;
            
            i = i+1;
        end
    end
    
    %% back to outer loop
    
    % selection of best parameter
    [params_sort, sort_idx] = sort(params_list(:, 4),'descend');
    params2 = params_list(:, 1:3);
    params2 = params2(sort_idx, :);
    params_list = [params2 params_sort];
    
    bestC = params_list(1,2);
    bestNu = params_list(1,3);
    
    options = ['-s 1 -t 2 -c ' num2str(bestC) ' -n ' num2str(bestNu) ' -q']; % Without cv; '-q' to disable screen output
    model = svmtrain(outcome_train, sparse(double(features_train)), options); % Need structure for prediction
    
    [results, predicted_label, accuracy, ~] = evalc('svmpredict(outcome_test, sparse(double(features_test)), model)');
    num_correct = extractBetween(results, '% (', ') (classification)');
    
    model_struc(i_oL).iteration_oL = i_oL;
    model_struc(i_oL).subjects_test_oL = subjects_test;
    model_struc(i_oL).bestC_iL = bestC;
    model_struc(i_oL).bestNu_iL = bestNu;
    model_struc(i_oL).outcome_test_oL = outcome_test;
    model_struc(i_oL).predicted_label_oL = predicted_label;
    model_struc(i_oL).group_sizes_train_oL = group_size_train;
    model_struc(i_oL).group_size_test_oL = group_size_test;
    model_struc(i_oL).params_list_iL = params_list;
    model_struc(i_oL).cv_iL = cv_struc;
    model_struc(i_oL).accuracy_oL = accuracy(1);
    
end

save('nestedCV_model_4class.mat','model_struc');

%% Average

for n = 1:10
    
    pred_all = [];
    test_all = [];
    subj_all = [];
    for s = 1:size(model_struc, 2)
        pred_all = [pred_all; model_struc(s).predicted_label_oL];
        test_all = [test_all; model_struc(s).outcome_test_oL];
        subj_all = [subj_all; model_struc(s).subjects_test_oL];
    end
    
    squRes = (test_all - pred_all).^2;
    squTotal = (test_all - mean(test_all)).^2;
    RSquare_all = 1 - (sum(squRes)/sum(squTotal));
    
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
    
end

pred_mean = mean([mat_sort.pred_sort], 2);
% squRes = (test_sort - pred_mean).^2;
% squTotal = (test_sort - mean(test_sort)).^2;
% RSquare_av = 1 - (sum(squRes)/sum(squTotal));

%%

acc = zeros(206,1);

for i = 1:206
    if isequal(pred_mean(i),test_sort(i))
        acc(i) = 1;
    end
end

acc_score = sum(acc(:))/206;

%%

pred_results = struct;
pred_results.mat_all = mat_all;
pred_results.mat_sort = mat_sort;
pred_results.acc_score = acc_score;

save('avgmodel_pred_results_4class.mat', 'pred_results');