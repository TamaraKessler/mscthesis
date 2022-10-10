
clear all; clc;
addpath('D:\MATLAB\libsvm-3.24\libsvm-3.24\matlab');
addpath('C:\Users\Sperber\Desktop\Sex Diff Scripts\OGs');
datadir = 'D:\Tamara\Diskonnektionen\LibSVM';
cd(datadir);

test_inst = readmatrix(fullfile(datadir, 'instance_matrix.csv'));
testsparse_inst = sparse(test_inst);

[N, M] = size(testsparse_inst);

test_lab = readmatrix(fullfile(datadir, 'labelsFM.csv'));

%%

defaultStream = RandStream('mlfg6331_64');
log2C = -5:1:14; % Ascending to pick the biggest value before reaching plateau
Nu = 0.01:0.05:1;
%model_version = 1;

model_struc = struct;

k = 5;

oL = k;
iL = length(Nu);
results = struct;

allpnum = 1:N;
allpnum = allpnum';

idx_data = (1:size(testsparse_inst, 1))';
idx_perm = randperm(defaultStream, length(idx_data))';
data_perm = testsparse_inst(idx_perm, :);
lab_perm = test_lab(idx_perm, :);
num_perm = allpnum(idx_perm, :);
    
split_numb = size(data_perm, 1) / k;
group_size_mean = floor(split_numb);
rest_numb = size(data_perm, 1) - (k(1) * group_size_mean);
group_sizes_v = zeros(k(1), 1);
group_sizes_v(1:(k(1)-rest_numb), 1) = group_size_mean;
group_sizes_v((k(1)-rest_numb+1):end, 1) = group_size_mean + 1;
if sum(group_sizes_v) ~= size(data_perm, 1) % Check
    fprintf('\nERROR. Wrong group sizes. \n');
    return
end

for i_oL = 1:k
    
    group_size_test = group_sizes_v(i_oL, 1);
    group_size_train = sum(group_sizes_v)-group_sizes_v(i_oL, 1);
    features_test = data_perm(1 : group_size_test, :); % 1 part serves as test set
    outcome_test = lab_perm(1 : group_size_test, :);
    subjects_test = num_perm(1 : group_size_test, :);
    features_train = data_perm(group_size_test+1 : end, :); % k-1 parts merged to one training set
    outcome_train = lab_perm(group_size_test+1 : end, :);
    subjects_train = num_perm(group_size_test+1 : end, :);
    
    % Move current test-set to end of matrix (preparation for next
    % iteration)
    data_perm = [data_perm(group_size_test+1 : end, :); ...
        data_perm(1:group_size_test, :)];
    lab_perm = [lab_perm(group_size_test+1 : end, :); ...
        lab_perm(1:group_size_test, :)];
    num_perm = [num_perm(group_size_test+1 : end, :); ...
        num_perm(1:group_size_test, :)];
    
    
    res = zeros(iL);
    
    for i_C = 1:iL
        
        for i_Nu = 1:iL
            
            if Nu(i_Nu)>=0.95
                options =  ['-s 1 -t 1 -c ' num2str(2^log2C(i_C)) ' -n ' num2str(Nu(i_Nu)-0.01) ' -b 1'];
            else
                options =  ['-s 1 -t 1 -c ' num2str(2^log2C(i_C)) ' -n ' num2str(Nu(i_Nu)) ' -b 1'];
            end
            model = svmtrain(outcome_train, features_train, options);
            [predicted_label, accuracy, ~] = svmpredict(outcome_test, features_test, model, '-b 1');
            
            res(i_C,i_Nu) = accuracy(1);
            
        end % i_Nu
        
    end %iL
    
    max_acc = max(res(1,:));
    idc_acc = find(res(1,:)==max_acc);
    best_nu = Nu(idc_acc);
    best_log2c = log2C(idc_acc);
    
%     if length(best_nu)>1
%         best_nu = best_nu(1);
%         best_log2c = best_log2c(1);
%     end

    model_struc(i_oL).iteration_oL = i_oL;
    model_struc(i_oL).subjects_test_oL = subjects_test;
    model_struc(i_oL).bestC_iL = best_log2c;
    model_struc(i_oL).bestNu_iL = best_nu;
    %model_struc(i_oL).bestMAE_iL = bestMAE;
    %model_struc(i_oL).model_MAE_oL = MAE_oL;
    model_struc(i_oL).outcome_test_oL = outcome_test;
    model_struc(i_oL).predicted_label_oL = predicted_label;
    model_struc(i_oL).max_accuracy = max_acc;
    
end

% pred_all = [];
% test_all = [];
% subj_all = [];
% for s = 1:size(model_struc, 2)
%     pred_all = [pred_all; model_struc(s).predicted_label_oL];
%     test_all = [test_all; model_struc(s).outcome_test_oL];
%     subj_all = [subj_all; model_struc(s).subjects_test_oL];
% end



%%

n = N;
l = n/2;
trainidc = randperm(n,l);

allnum = 1:N;
testidc = setdiff(allnum,trainidc);

trainIndex = zeros(N,1); trainIndex(trainidc) = 1;
testIndex = zeros(N,1); testIndex(testidc) = 1;
trainData = testsparse_inst(trainIndex==1,:);
trainLabel = test_lab(trainIndex==1,:);
testData = testsparse_inst(testIndex==1,:);
testLabel = test_lab(testIndex==1,:);

options = ['-s 1 -t 1 -c 5 -n 0.01 -b 1']

model = svmtrain(trainLabel, trainData, options);
[predict_label, accuracy, prob_values] = svmpredict(testLabel, testData, model, '-b 1');

% [predicted_label, accuracy, decision_values/prob_estimates] = svmpredict(testing_label_vector, testing_instance_matrix, model [, 'libsvm_options']);

