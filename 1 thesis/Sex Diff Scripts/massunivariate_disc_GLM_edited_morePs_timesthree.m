%%
% A L L
%% White Matter Disconnection Analysis for disconnection-symptom mapping
% via GLM for continuous data

% by Christoph Sperber

% This script requires ROI-to-ROI disconnection matrices of patients (e.g.
% obtained from NeMo, BCB or LQT) and behavioural scores. For each
% ROI-to-ROI connection, a GLM is computed to test the relation between
% disconnection and behaviour. This is a classical mass-univariate approach
% to test for brain-behaviour relations. Statistical significance is
% assessed by maximum statistic permutation - this approach implements an
% exact FWER correction to control for multiple comparisons

% many ROI-to-ROI connections are
% likely irrelavant because they are physiologically non-existent or never affected
% in the lesion sample, and the data are (almost)
% always zero. Therefore, before the actual analysis the
% data are masked to remove any ROI-to-ROI connections that are rarely affected.
% This cutoff is arbitrary, just like the minimum-lesion-affection
% criterion in VLSM ('we tested only voxels affected by at least n
% lesions'), but choosing a cutoff a priori (!) is
% fine, just like in VLSM (see Sperber & Karnath, 2017 in HBM)

clear
clc
start_time = clock;
f = waitbar(0,'Loading Data'); % this initialises a bar that depicts the progress of permutations

% define number of permutations in max stat permutation. Importantly, choose a large number (at least 1.000, better 10.000 or even
% more), and make sure that the final p-threshold can be assessed with the
% resolution of p-values provided by the permutation (e.g. mapping results
% at p<0.0001 while having only 1000 permutations is non-sense, as the
% lowest possible p-values are 0.001 and 0)
perms =50000;

% choose the 1-p statistical threshold (p < 0.05 would be 0.95)
stat_threshold = [0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999];

% read behavioural data
behaviour =  xlsread('D:\Tamara\DisconnectionMaps\Parcel_Disconnection\all\renamed\all_behavioural_discon_neg.xlsx');
% this line just reads the number from the xls. If future matlab version
% should change anything with the xlsread or if you have any issues: just
% make sure that 'behaviour' is a n-by-1 vector that contains the
% behavioural scores
% IMPORTANT: check the polarity of the data - in the current version, the algorithms assume
% higher scores mean BETTER (=less pathological) performance , therefore the final statistics are
% +/- inversed to allow easier interpretation. Thus, damage to areas with positive
% statistics is associated with more deficit
% IMPORTANT: to simplify the scripting, it is assumed that the
% behavioural scores are in the same order as the disconnection files

% Set default stream for pseudorandomisation
defaultStream = RandStream('mlfg6331_64');

%% read imaging data
folder = 'D:\Tamara\DisconnectionMaps\Parcel_Disconnection\all\renamed'; % this folder contains all relevant disconnection .mat files and NO other .mat files
fileList = struct2cell(dir(fullfile(folder, '*.mat')));
% short check: mat-files need to be in the same order as in the
% behavioural file; ideally, double check this here!
temp=load(strcat(fileList{2,1}, '\', fileList{1,1})).pct_sdc_matrix; %load first file to get matrix size
% IMPORTANT: depending on the format of the data you have to modify this
% line. This line works with the LQT output, which is a struct that
% contains the a variable 'pct_sdc_matrix'. The '.' is the syntax to access
% an element within a struct. If you open a text file, an excel file or a
% matlab file without such a struct, you have to adapt this line

images_2d = zeros(length(temp),length(temp),length(fileList)); %create blank for 2d images

for i=1:length(fileList)
    images_2d(:,:,i)=load(strcat(fileList{2,i}, '\', fileList{1,i})).pct_sdc_matrix;
    % IMPORTANT: same thing as in the previous comment!
end


%% create a 2d mask for the imaging data
mask=zeros(length(temp));
% the matrix is a symmetric matrix, so first remove all elements in the
% diagonal and below
for i=1:length(temp)
    for j=1:length(temp)
        if i<j
            mask(i,j)=1;
        end
    end
end

% then we need to remove rarely/never affected connections
images_2d_binary=double(images_2d>0);% binarise the image, every disconnection >0 is set to 1
image_2d_sum=sum(images_2d_binary,3);
image_2d_sum=image_2d_sum.*mask;

% get an overview on the number of patients with a disconnection for each
% ROI to ROI connection
% temp=sort(image_2d_sum(:),'descend');
% in 'temp', you can check how many connections are affected in at least x patients

% you can choose to only test connections that are affected in at least threshold = x
% patients (I set it to 5, but feel free to change it)
threshold = 55;
mask(image_2d_sum<threshold)=0;
mask_vect=mask(:);

% vectorise the 2d images
images_vect=zeros(length(fileList),sum(mask(:))); %create blank
for i=1:length(fileList)
    temp=images_2d(:,:,i);
    temp=temp(:);
    images_vect(i,:)=temp(mask_vect==1);
end

%% do the GLM with the original data
waitbar(0,f,'Starting GLM');
stat_vect_orig=zeros(sum(mask(:)),1);
p_vect_orig=zeros(sum(mask(:)),1);

for i=1:length(stat_vect_orig)
    [b,dev,stats]=glmfit(images_vect(:,i),behaviour);
    stat_vect_orig(i)=-(stats.t(2)); % the minus sign changes the direction, see comment above on polarity
    p_vect_orig(i)=stats.p(2);
end

%% and now the same GLM for permuted data. In each run, only the maximum
% statistic is saved
permutation_max_stats=zeros(perms,1); % create a blank
stat_vect_perm_temp=zeros(sum(mask(:)),1);

for j=1:perms
    waitbar(j/perms,f,'Doing GLM permutations');

    % Pseudorandomisation
    defaultStream.Substream = j;
    rand_behaviour=behaviour(randperm(defaultStream,length(behaviour)));

    % Randomisation
    %rand_behaviour=behaviour(randperm(length(behaviour))); % permute behavioural data

    for i=1:length(stat_vect_orig)
        [b,dev,stats]=glmfit(images_vect(:,i),rand_behaviour);
        stat_vect_perm_temp(i)=-(stats.t(2)); % the minus sign changes the direction, see comment above on polarity
    end
    permutation_max_stats(j)=max(stat_vect_perm_temp);
end

% the permutation-derived max statistics are now ordered and the relevent
% percentile is identified (e.g. the 95th percentile with a p value of
% 0.05)
permutation_max_stats = sort(permutation_max_stats);
% this is the statistical threshold. Every imaging feature in the original
% statistical map larger than this threshold is significant at the chosen
% p-level

for i_thresh = 1:length(stat_threshold)

    permutation_threshold = permutation_max_stats(round(perms*stat_threshold(i_thresh)));

    % additional FDR analysis, using the Benjamini Yekuteli method via fdr_bh.m
    % (make sure to add the function to your matlab path!)
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_vect_orig,1-stat_threshold(i_thresh),'dep');
    % move the result back into original matrix space
    results_2d_fdr = zeros(length(mask));
    temp=mask(:);
    temp(temp==1)=h;
    results_2d_fdr = reshape(temp,[length(mask),length(mask)]);


    % move the result back into original matrix space
    results_2d = zeros(length(mask));

    temp=mask(:);
    temp(temp==1)=p_vect_orig;
    p_values_orig = reshape(temp,[length(mask),length(mask)]); % in this matrix, every tested voxel contains its original glm p-value

    temp=mask(:);
    temp(temp==1)=stat_vect_orig;
    results_2d = reshape(temp,[length(mask),length(mask)]); % in this matrix, every tested voxel contains its original glm statistic
    results_2d_thresholded = double(results_2d>permutation_threshold); % in this matrix, every significant voxel at the chosen p is 1, else 0

    % a small hint to make this matrix more user-friendly: a large symmetric matrix can contain a lot
    % of cells, you wont get a good overview. Instead, transform it to
    % sparse format and get a full list of all significant cells
    results_thresholded_view = sparse(results_2d_thresholded);
    % and now just type:
    % results_thresholded_view
    % and see all matrix cells with significant results.


    % create a results file and add all kinds of informative stuff
    results = struct;
    results.perm_num = perms; % permutation number
    results.p_level = 1-stat_threshold(i_thresh); % user-defined p-value
    results.input_data = fileList; % files used in the analysis
    results.mask_tested_disconns = mask; % these connections were included in the analysis
    results.max_stat_threshold = permutation_threshold; % this is the FWER-corrected t-statistic threshold for significance
    results.glm_results_vect = stat_vect_orig; % these are the feature-wise t-statistics of the original (uncorrected) analysis
    results.time_finished = clock;
    results.perm_max_stats = permutation_max_stats; % this is a vector of all maximum statistics for all permutations
    results.orig_p_values=p_values_orig; % these are the p-values in the original analysis
    results.resulting_glm_t_statistic = results_2d; % these are the feature-wise statistics of the original (uncorrected) analysis in original 2d format
    results.resulting_statistics_binary = results_2d_thresholded; % these are all matrix cells that are significant after fwer perumtation correction (binary, 1=significant)
    results.fdr_threshold = crit_p; % this is the critical p for fdr control
    results.fdr_controlled_binary=results_2d_fdr; % this binarily shows all matrix cells that are siginficant after FDR control (binary, 1=significant)
    results.randSeed = defaultStream;

    % to simplify everything: if you want
    % a) permutation-wise fwer by maximum statistics permuation
    % --> either use 'results.resulting_statistics_binary', which shows
    % binarily which matrix cells are significant, or use
    % 'results.resulting_glm_t_statistic' with a threshold defined by
    % 'results.max_stat_threshold' if you want to visualise/use the range of
    % statistical values above the FWER corrected threshold
    % b) FDR correction
    % --> use 'results.orig_p_values' and apply the threshold in
    % 'results.fdr_threshold', or use 'results.fdr_controlled_binary' for the
    % binary map of all matrix cells that are significant

    saveloc  = 'D:\Tamara\DisconnectionMaps\Parcel_Disconnection\results_GLM\all\thresh25perc';
    savename = strcat('results_disconn_map_GLM_p_',num2str(10000-(stat_threshold(i_thresh)*10000)), '_', ...
        num2str(start_time(1)), '_' , num2str(start_time(2)), '_' , num2str(start_time(3)), '_' , num2str(start_time(4)), '_' , num2str(start_time(5)));
    full_savename = fullfile(saveloc, savename);
    save(full_savename,'results');

    dlmwrite([full_savename '.edge'], results.resulting_statistics_binary, 'delimiter', '\t', 'precision', 4);

end

% close the waitbar
close(f)
fprintf('Done')

%%
% F E M A L E 
%%

%% White Matter Disconnection Analysis for disconnection-symptom mapping
% via GLM for continuous data

% by Christoph Sperber

% This script requires ROI-to-ROI disconnection matrices of patients (e.g.
% obtained from NeMo, BCB or LQT) and behavioural scores. For each
% ROI-to-ROI connection, a GLM is computed to test the relation between
% disconnection and behaviour. This is a classical mass-univariate approach
% to test for brain-behaviour relations. Statistical significance is
% assessed by maximum statistic permutation - this approach implements an
% exact FWER correction to control for multiple comparisons

% many ROI-to-ROI connections are
% likely irrelavant because they are physiologically non-existent or never affected
% in the lesion sample, and the data are (almost)
% always zero. Therefore, before the actual analysis the
% data are masked to remove any ROI-to-ROI connections that are rarely affected.
% This cutoff is arbitrary, just like the minimum-lesion-affection
% criterion in VLSM ('we tested only voxels affected by at least n
% lesions'), but choosing a cutoff a priori (!) is
% fine, just like in VLSM (see Sperber & Karnath, 2017 in HBM)

clear
clc
start_time = clock;
f = waitbar(0,'Loading Data'); % this initialises a bar that depicts the progress of permutations

% define number of permutations in max stat permutation. Importantly, choose a large number (at least 1.000, better 10.000 or even
% more), and make sure that the final p-threshold can be assessed with the
% resolution of p-values provided by the permutation (e.g. mapping results
% at p<0.0001 while having only 1000 permutations is non-sense, as the
% lowest possible p-values are 0.001 and 0)
perms =50000;

% choose the 1-p statistical threshold (p < 0.05 would be 0.95)
stat_threshold = [0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999];

% read behavioural data
behaviour =  xlsread('D:\Tamara\DisconnectionMaps\Parcel_Disconnection\all_female\all_female_behavioural_discon_neg.xlsx');
% this line just reads the number from the xls. If future matlab version
% should change anything with the xlsread or if you have any issues: just
% make sure that 'behaviour' is a n-by-1 vector that contains the
% behavioural scores
% IMPORTANT: check the polarity of the data - in the current version, the algorithms assume
% higher scores mean BETTER (=less pathological) performance , therefore the final statistics are
% +/- inversed to allow easier interpretation. Thus, damage to areas with positive
% statistics is associated with more deficit
% IMPORTANT: to simplify the scripting, it is assumed that the
% behavioural scores are in the same order as the disconnection files

% Set default stream for pseudorandomisation
defaultStream = RandStream('mlfg6331_64');

%% read imaging data
folder = 'D:\Tamara\DisconnectionMaps\Parcel_Disconnection\all_female'; % this folder contains all relevant disconnection .mat files and NO other .mat files
fileList = struct2cell(dir(fullfile(folder, '*.mat')));
% short check: mat-files need to be in the same order as in the
% behavioural file; ideally, double check this here!
temp=load(strcat(fileList{2,1}, '\', fileList{1,1})).pct_sdc_matrix; %load first file to get matrix size
% IMPORTANT: depending on the format of the data you have to modify this
% line. This line works with the LQT output, which is a struct that
% contains the a variable 'pct_sdc_matrix'. The '.' is the syntax to access
% an element within a struct. If you open a text file, an excel file or a
% matlab file without such a struct, you have to adapt this line

images_2d = zeros(length(temp),length(temp),length(fileList)); %create blank for 2d images

for i=1:length(fileList)
    images_2d(:,:,i)=load(strcat(fileList{2,i}, '\', fileList{1,i})).pct_sdc_matrix;
    % IMPORTANT: same thing as in the previous comment!
end


%% create a 2d mask for the imaging data
mask=zeros(length(temp));
% the matrix is a symmetric matrix, so first remove all elements in the
% diagonal and below
for i=1:length(temp)
    for j=1:length(temp)
        if i<j
            mask(i,j)=1;
        end
    end
end

% then we need to remove rarely/never affected connections
images_2d_binary=double(images_2d>0);% binarise the image, every disconnection >0 is set to 1
image_2d_sum=sum(images_2d_binary,3);
image_2d_sum=image_2d_sum.*mask;

% get an overview on the number of patients with a disconnection for each
% ROI to ROI connection
% temp=sort(image_2d_sum(:),'descend');
% in 'temp', you can check how many connections are affected in at least x patients

% you can choose to only test connections that are affected in at least threshold = x
% patients (I set it to 5, but feel free to change it)
threshold = 25;
mask(image_2d_sum<threshold)=0;
mask_vect=mask(:);

% vectorise the 2d images
images_vect=zeros(length(fileList),sum(mask(:))); %create blank
for i=1:length(fileList)
    temp=images_2d(:,:,i);
    temp=temp(:);
    images_vect(i,:)=temp(mask_vect==1);
end

%% do the GLM with the original data
waitbar(0,f,'Starting GLM');
stat_vect_orig=zeros(sum(mask(:)),1);
p_vect_orig=zeros(sum(mask(:)),1);

for i=1:length(stat_vect_orig)
    [b,dev,stats]=glmfit(images_vect(:,i),behaviour);
    stat_vect_orig(i)=-(stats.t(2)); % the minus sign changes the direction, see comment above on polarity
    p_vect_orig(i)=stats.p(2);
end

%% and now the same GLM for permuted data. In each run, only the maximum
% statistic is saved
permutation_max_stats=zeros(perms,1); % create a blank
stat_vect_perm_temp=zeros(sum(mask(:)),1);

for j=1:perms
    waitbar(j/perms,f,'Doing GLM permutations');

    % Pseudorandomisation
    defaultStream.Substream = j;
    rand_behaviour=behaviour(randperm(defaultStream,length(behaviour)));

    % Randomisation
    %rand_behaviour=behaviour(randperm(length(behaviour))); % permute behavioural data

    for i=1:length(stat_vect_orig)
        [b,dev,stats]=glmfit(images_vect(:,i),rand_behaviour);
        stat_vect_perm_temp(i)=-(stats.t(2)); % the minus sign changes the direction, see comment above on polarity
    end
    permutation_max_stats(j)=max(stat_vect_perm_temp);
end

% the permutation-derived max statistics are now ordered and the relevent
% percentile is identified (e.g. the 95th percentile with a p value of
% 0.05)
permutation_max_stats = sort(permutation_max_stats);
% this is the statistical threshold. Every imaging feature in the original
% statistical map larger than this threshold is significant at the chosen
% p-level

for i_thresh = 1:length(stat_threshold)

    permutation_threshold = permutation_max_stats(round(perms*stat_threshold(i_thresh)));

    % additional FDR analysis, using the Benjamini Yekuteli method via fdr_bh.m
    % (make sure to add the function to your matlab path!)
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_vect_orig,1-stat_threshold(i_thresh),'dep');
    % move the result back into original matrix space
    results_2d_fdr = zeros(length(mask));
    temp=mask(:);
    temp(temp==1)=h;
    results_2d_fdr = reshape(temp,[length(mask),length(mask)]);


    % move the result back into original matrix space
    results_2d = zeros(length(mask));

    temp=mask(:);
    temp(temp==1)=p_vect_orig;
    p_values_orig = reshape(temp,[length(mask),length(mask)]); % in this matrix, every tested voxel contains its original glm p-value

    temp=mask(:);
    temp(temp==1)=stat_vect_orig;
    results_2d = reshape(temp,[length(mask),length(mask)]); % in this matrix, every tested voxel contains its original glm statistic
    results_2d_thresholded = double(results_2d>permutation_threshold); % in this matrix, every significant voxel at the chosen p is 1, else 0

    % a small hint to make this matrix more user-friendly: a large symmetric matrix can contain a lot
    % of cells, you wont get a good overview. Instead, transform it to
    % sparse format and get a full list of all significant cells
    results_thresholded_view = sparse(results_2d_thresholded);
    % and now just type:
    % results_thresholded_view
    % and see all matrix cells with significant results.


    % create a results file and add all kinds of informative stuff
    results = struct;
    results.perm_num = perms; % permutation number
    results.p_level = 1-stat_threshold(i_thresh); % user-defined p-value
    results.input_data = fileList; % files used in the analysis
    results.mask_tested_disconns = mask; % these connections were included in the analysis
    results.max_stat_threshold = permutation_threshold; % this is the FWER-corrected t-statistic threshold for significance
    results.glm_results_vect = stat_vect_orig; % these are the feature-wise t-statistics of the original (uncorrected) analysis
    results.time_finished = clock;
    results.perm_max_stats = permutation_max_stats; % this is a vector of all maximum statistics for all permutations
    results.orig_p_values=p_values_orig; % these are the p-values in the original analysis
    results.resulting_glm_t_statistic = results_2d; % these are the feature-wise statistics of the original (uncorrected) analysis in original 2d format
    results.resulting_statistics_binary = results_2d_thresholded; % these are all matrix cells that are significant after fwer perumtation correction (binary, 1=significant)
    results.fdr_threshold = crit_p; % this is the critical p for fdr control
    results.fdr_controlled_binary=results_2d_fdr; % this binarily shows all matrix cells that are siginficant after FDR control (binary, 1=significant)
    results.randSeed = defaultStream;

    % to simplify everything: if you want
    % a) permutation-wise fwer by maximum statistics permuation
    % --> either use 'results.resulting_statistics_binary', which shows
    % binarily which matrix cells are significant, or use
    % 'results.resulting_glm_t_statistic' with a threshold defined by
    % 'results.max_stat_threshold' if you want to visualise/use the range of
    % statistical values above the FWER corrected threshold
    % b) FDR correction
    % --> use 'results.orig_p_values' and apply the threshold in
    % 'results.fdr_threshold', or use 'results.fdr_controlled_binary' for the
    % binary map of all matrix cells that are significant

    saveloc  = 'D:\Tamara\DisconnectionMaps\Parcel_Disconnection\results_GLM\female\thresh25perc';
    savename = strcat('results_disconn_map_GLM_p_',num2str(10000-(stat_threshold(i_thresh)*10000)), '_', ...
        num2str(start_time(1)), '_' , num2str(start_time(2)), '_' , num2str(start_time(3)), '_' , num2str(start_time(4)), '_' , num2str(start_time(5)));
    full_savename = fullfile(saveloc, savename);
    save(full_savename,'results');

    dlmwrite([full_savename '.edge'], results.resulting_statistics_binary, 'delimiter', '\t', 'precision', 4);

end

% close the waitbar
close(f)
fprintf('Done')

%%
% M A L E
%%
%% White Matter Disconnection Analysis for disconnection-symptom mapping
% via GLM for continuous data

% by Christoph Sperber

% This script requires ROI-to-ROI disconnection matrices of patients (e.g.
% obtained from NeMo, BCB or LQT) and behavioural scores. For each
% ROI-to-ROI connection, a GLM is computed to test the relation between
% disconnection and behaviour. This is a classical mass-univariate approach
% to test for brain-behaviour relations. Statistical significance is
% assessed by maximum statistic permutation - this approach implements an
% exact FWER correction to control for multiple comparisons

% many ROI-to-ROI connections are
% likely irrelavant because they are physiologically non-existent or never affected
% in the lesion sample, and the data are (almost) 
% always zero. Therefore, before the actual analysis the
% data are masked to remove any ROI-to-ROI connections that are rarely affected. 
% This cutoff is arbitrary, just like the minimum-lesion-affection
% criterion in VLSM ('we tested only voxels affected by at least n
% lesions'), but choosing a cutoff a priori (!) is
% fine, just like in VLSM (see Sperber & Karnath, 2017 in HBM)

clear
clc
start_time = clock;
f = waitbar(0,'Loading Data'); % this initialises a bar that depicts the progress of permutations

% define number of permutations in max stat permutation. Importantly, choose a large number (at least 1.000, better 10.000 or even
% more), and make sure that the final p-threshold can be assessed with the
% resolution of p-values provided by the permutation (e.g. mapping results
% at p<0.0001 while having only 1000 permutations is non-sense, as the
% lowest possible p-values are 0.001 and 0)
perms =50000;

% choose the 1-p statistical threshold (p < 0.05 would be 0.95)
stat_threshold = [0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999];

% read behavioural data
behaviour =  xlsread('D:\Tamara\DisconnectionMaps\Parcel_Disconnection\all_male\all_male_behavioural_neg.xlsx');
% this line just reads the number from the xls. If future matlab version
% should change anything with the xlsread or if you have any issues: just
% make sure that 'behaviour' is a n-by-1 vector that contains the
% behavioural scores
% IMPORTANT: check the polarity of the data - in the current version, the algorithms assume
% higher scores mean BETTER (=less pathological) performance , therefore the final statistics are
% +/- inversed to allow easier interpretation. Thus, damage to areas with positive
% statistics is associated with more deficit
% IMPORTANT: to simplify the scripting, it is assumed that the
% behavioural scores are in the same order as the disconnection files

% Set default stream for pseudorandomisation
defaultStream = RandStream('mlfg6331_64');

%% read imaging data
folder = 'D:\Tamara\DisconnectionMaps\Parcel_Disconnection\all_male'; % this folder contains all relevant disconnection .mat files and NO other .mat files
fileList = struct2cell(dir(fullfile(folder, '*.mat')));
% short check: mat-files need to be in the same order as in the
% behavioural file; ideally, double check this here!
temp=load(strcat(fileList{2,1}, '\', fileList{1,1})).pct_sdc_matrix; %load first file to get matrix size
% IMPORTANT: depending on the format of the data you have to modify this
% line. This line works with the LQT output, which is a struct that
% contains the a variable 'pct_sdc_matrix'. The '.' is the syntax to access
% an element within a struct. If you open a text file, an excel file or a
% matlab file without such a struct, you have to adapt this line

images_2d = zeros(length(temp),length(temp),length(fileList)); %create blank for 2d images

for i=1:length(fileList)
    images_2d(:,:,i)=load(strcat(fileList{2,i}, '\', fileList{1,i})).pct_sdc_matrix;
    % IMPORTANT: same thing as in the previous comment!
end


%% create a 2d mask for the imaging data
mask=zeros(length(temp));
% the matrix is a symmetric matrix, so first remove all elements in the
% diagonal and below
for i=1:length(temp)
    for j=1:length(temp)
        if i<j
            mask(i,j)=1;
        end
    end
end

% then we need to remove rarely/never affected connections
images_2d_binary=double(images_2d>0);% binarise the image, every disconnection >0 is set to 1
image_2d_sum=sum(images_2d_binary,3);
image_2d_sum=image_2d_sum.*mask;

% get an overview on the number of patients with a disconnection for each
% ROI to ROI connection
% temp=sort(image_2d_sum(:),'descend');
% in 'temp', you can check how many connections are affected in at least x patients

% you can choose to only test connections that are affected in at least threshold = x
% patients (I set it to 5, but feel free to change it)
threshold = 30;
mask(image_2d_sum<threshold)=0;
mask_vect=mask(:);

% vectorise the 2d images
images_vect=zeros(length(fileList),sum(mask(:))); %create blank
for i=1:length(fileList)
    temp=images_2d(:,:,i);
    temp=temp(:);
    images_vect(i,:)=temp(mask_vect==1);
end

%% do the GLM with the original data
waitbar(0,f,'Starting GLM');
stat_vect_orig=zeros(sum(mask(:)),1);
p_vect_orig=zeros(sum(mask(:)),1);

for i=1:length(stat_vect_orig)
    [b,dev,stats]=glmfit(images_vect(:,i),behaviour);
    stat_vect_orig(i)=-(stats.t(2)); % the minus sign changes the direction, see comment above on polarity
    p_vect_orig(i)=stats.p(2);
end

%% and now the same GLM for permuted data. In each run, only the maximum
% statistic is saved
permutation_max_stats=zeros(perms,1); % create a blank
stat_vect_perm_temp=zeros(sum(mask(:)),1);

for j=1:perms
    waitbar(j/perms,f,'Doing GLM permutations');
    
    % Pseudorandomisation
    defaultStream.Substream = j;
    rand_behaviour=behaviour(randperm(defaultStream,length(behaviour)));
    
    % Randomisation
    %rand_behaviour=behaviour(randperm(length(behaviour))); % permute behavioural data
    
    for i=1:length(stat_vect_orig)
        [b,dev,stats]=glmfit(images_vect(:,i),rand_behaviour);
        stat_vect_perm_temp(i)=-(stats.t(2)); % the minus sign changes the direction, see comment above on polarity
    end
    permutation_max_stats(j)=max(stat_vect_perm_temp);
end

% the permutation-derived max statistics are now ordered and the relevent
% percentile is identified (e.g. the 95th percentile with a p value of
% 0.05)
permutation_max_stats = sort(permutation_max_stats);
% this is the statistical threshold. Every imaging feature in the original
% statistical map larger than this threshold is significant at the chosen
% p-level

for i_thresh = 1:length(stat_threshold)

    permutation_threshold = permutation_max_stats(round(perms*stat_threshold(i_thresh)));

    % additional FDR analysis, using the Benjamini Yekuteli method via fdr_bh.m
    % (make sure to add the function to your matlab path!)
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p_vect_orig,1-stat_threshold(i_thresh),'dep');
    % move the result back into original matrix space
    results_2d_fdr = zeros(length(mask));
    temp=mask(:);
    temp(temp==1)=h;
    results_2d_fdr = reshape(temp,[length(mask),length(mask)]);
    
    
    % move the result back into original matrix space
    results_2d = zeros(length(mask));
    
    temp=mask(:);
    temp(temp==1)=p_vect_orig;
    p_values_orig = reshape(temp,[length(mask),length(mask)]); % in this matrix, every tested voxel contains its original glm p-value
    
    temp=mask(:);
    temp(temp==1)=stat_vect_orig;
    results_2d = reshape(temp,[length(mask),length(mask)]); % in this matrix, every tested voxel contains its original glm statistic
    results_2d_thresholded = double(results_2d>permutation_threshold); % in this matrix, every significant voxel at the chosen p is 1, else 0
    
    % a small hint to make this matrix more user-friendly: a large symmetric matrix can contain a lot
    % of cells, you wont get a good overview. Instead, transform it to
    % sparse format and get a full list of all significant cells
    results_thresholded_view = sparse(results_2d_thresholded);
    % and now just type:
    % results_thresholded_view
    % and see all matrix cells with significant results.
    
    
    % create a results file and add all kinds of informative stuff
    results = struct;
    results.perm_num = perms; % permutation number
    results.p_level = 1-stat_threshold(i_thresh); % user-defined p-value
    results.input_data = fileList; % files used in the analysis
    results.mask_tested_disconns = mask; % these connections were included in the analysis
    results.max_stat_threshold = permutation_threshold; % this is the FWER-corrected t-statistic threshold for significance
    results.glm_results_vect = stat_vect_orig; % these are the feature-wise t-statistics of the original (uncorrected) analysis
    results.time_finished = clock;
    results.perm_max_stats = permutation_max_stats; % this is a vector of all maximum statistics for all permutations
    results.orig_p_values=p_values_orig; % these are the p-values in the original analysis
    results.resulting_glm_t_statistic = results_2d; % these are the feature-wise statistics of the original (uncorrected) analysis in original 2d format
    results.resulting_statistics_binary = results_2d_thresholded; % these are all matrix cells that are significant after fwer perumtation correction (binary, 1=significant)
    results.fdr_threshold = crit_p; % this is the critical p for fdr control
    results.fdr_controlled_binary=results_2d_fdr; % this binarily shows all matrix cells that are siginficant after FDR control (binary, 1=significant)
    results.randSeed = defaultStream;
    
    % to simplify everything: if you want
    % a) permutation-wise fwer by maximum statistics permuation
    % --> either use 'results.resulting_statistics_binary', which shows
    % binarily which matrix cells are significant, or use
    % 'results.resulting_glm_t_statistic' with a threshold defined by
    % 'results.max_stat_threshold' if you want to visualise/use the range of
    % statistical values above the FWER corrected threshold
    % b) FDR correction
    % --> use 'results.orig_p_values' and apply the threshold in
    % 'results.fdr_threshold', or use 'results.fdr_controlled_binary' for the
    % binary map of all matrix cells that are significant
    
    saveloc  = 'D:\Tamara\DisconnectionMaps\Parcel_Disconnection\results_GLM\male\thresh25perc';
    savename = strcat('results_disconn_map_GLM_p_',num2str(10000-(stat_threshold(i_thresh)*10000)), '_', ...
        num2str(start_time(1)), '_' , num2str(start_time(2)), '_' , num2str(start_time(3)), '_' , num2str(start_time(4)), '_' , num2str(start_time(5)));
    full_savename = fullfile(saveloc, savename);
    save(full_savename,'results');
    
    dlmwrite([full_savename '.edge'], results.resulting_statistics_binary, 'delimiter', '\t', 'precision', 4);

end

% close the waitbar
close(f)
fprintf('Done')
