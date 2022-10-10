%% White Matter Disconnection Analysis for number cognition deficits
% by Christoph Sperber for Stefan Smaczny
% This script requires ROI-to-ROI disconnection matrices of patients (e.g.
% obtained from NeMo or LQT) and their behavioural scores. For each
% ROI-to-ROI connection, a GLM is computed to test the relation between
% disconnection and behaviour. This is a classical mass-univariate approach
% to test for brain-behaviour relations. Statistical significance is
% assessed by maximum statistic permutation - this approach implements an
% exact FWER correction to control for multiple comparisons

% some important comments on the analysis strategy:
% First, many ROI-to-ROI connections are
% likely irrelavant because they are physiologically non-existent or never affected
% in the lesion sample, and the data are (almost) 
% always zero. Therefore, before the actual analysis the
% data are masked to remove any ROI-to-ROI connections that are rarely affected. 
% This cutoff is arbitrary, just like the minimum-lesion-affection
% criterion in VLSM ('we tested only voxels affected by at least n
% lesions'), but choosing a cutoff a priori (!) is
% fine, just like in VLSM
% Second, there might be a large amount of significant results (or might
% not, I don't know a priori). In VLSM, we
% sometimes look at 500.000 voxels, and 35.000 turn out to be significant.
% These voxels can then be interpreted by reference to a brain atlas.
% If in the current analysis we look at 500 disconnections and find 35 to
% be significant, this becomes an overly complex result, which will not
% allow any straightforward interpretation. Also, as always in lesion-brain
% inference, there will likely be a bunch of non-causal lesion-brain 
% associations involved (i.e. collateral damage; Sperber, 2020 in Cortex). Therefore, we start with a
% standard p-cutoff (0.05), but leave open the possibility to post-hoc
% adjust the p-value to a more conservative threshold. In this case, all
% results should be transparently reported with full p-values (e.g. in
% online/supplementary materials), but the final results & discussion will
% focus on a smaller set of results with a slightly more conservative p-cutoff that allows meaningful theoretical
% discussion; two further options are to either look at some specific
% connections/areas in the light of hypotheses, or to look if there are
% many disconnections to a single area, i.e. if there is some hub-like
% structure involved.
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
stat_threshold = 0.92;

% read behavioural data
behaviour =  xlsread('C:\Users\ssmaczny\Desktop\Lesion_Quantification_Toolkit\Zeug fuer Chris\Behavioural.xlsx');
% IMPORTANT: check the polarity of the data - in the multiplication task,
% higher scores mean BETTER performance, therefore the final statistics are
% +/- inversed to allow easier interpretation. Thus, damage to areas with positive
% statistics is associated with more deficit

%% read imaging data
folder = 'C:\Users\ssmaczny\Desktop\Lesion_Quantification_Toolkit\Zeug fuer Chris\Parcel Disconnection'; % this folder contains all disconnection mats and NO other mat files
fileList = struct2cell(dir(fullfile(folder, '*.mat')));
% short check: mat-files need to be in the same order as in the
% behavioural file; I did so and it looked fine
temp=load(strcat(fileList{2,1}, '\', fileList{1,1})).pct_sdc_matrix; %load first file to get matrix size
images_2d = zeros(length(temp),length(temp),length(fileList)); %create blank for 2d images

for i=1:length(fileList)
    images_2d(:,:,i)=load(strcat(fileList{2,i}, '\', fileList{1,i})).pct_sdc_matrix;
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
% the maximum number of patients with disconnection is 39; there are 2032
% connections affected in at least 10 patients and 955
% connections affected in at least 15 patients and

% I chose to only test connections that are affected in at least 15
% patients; feel free to change this to another value before you analyse
% the data
mask(image_2d_sum<15)=0;
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

for i=1:length(stat_vect_orig)
    [b,dev,stats]=glmfit(images_vect(:,i),behaviour);
    stat_vect_orig(i)=-(stats.t(2)); % the minus sign changes the direction, see comment above on polarity
end

%% and now the same GLM for permuted data. In each run, only the maximum
% statistic is saved
permutation_max_stats=zeros(perms,1); % create a blank
stat_vect_perm_temp=zeros(sum(mask(:)),1);

for j=1:perms
    waitbar(j/perms,f,'Doing GLM permutations');
    rand_behaviour=behaviour(randperm(length(behaviour))); % permute behavioural data
    for i=1:length(stat_vect_orig)
        [b,dev,stats]=glmfit(images_vect(:,i),rand_behaviour);
        stat_vect_perm_temp(i)=-(stats.t(2)); % the minus sign changes the direction, see comment above on polarity
    end
    permutation_max_stats(j)=max(stat_vect_perm_temp);
end

% the permutation-derived max statistics are now ordered and the 95th
% percentile is identified
permutation_max_stats = sort(permutation_max_stats);
% this is the statistical threshold. Every imaging feature in the original
% statistical map larger than this threshold is significant at the chosen
% p-level
permutation_threshold = permutation_max_stats(round(perms*stat_threshold));

% move the result back into original matrix space
results_2d = zeros(length(mask));

temp=mask(:);
temp(temp==1)=stat_vect_orig;
results_2d = reshape(temp,[length(mask),length(mask)]); % in this matrix, every tested voxel contains its original glm statistic
results_2d_thresholded = double(results_2d>permutation_threshold); % in this matrix, every significant voxel at the chosen p is 1, else 0

% a small hint to make this matrix more user-friendly: 246x246 cells are
% quite a lot, you wont get a good overview. Instead, transform it to
% sparse format and get a full list of all significant cells
results_thresholded_view = sparse(results_2d_thresholded);
% and now just type:
% results_thresholded_view
% and see all matrix cells with significant results. 


% create a results file and add all kinds of informative stuff
results = struct;
results.perm_num = perms;
results.p_level = 1-stat_threshold;
results.input_data = fileList;
results.mask_tested_disconns = mask;
results.max_stat_threshold = permutation_threshold;
results.glm_results_vect = stat_vect_orig;
results.time_finished = clock;
results.perm_max_stats = permutation_max_stats;
results.resulting_statistics = results_2d;
results.resulting_statistics_binary = results_2d_thresholded;

savename = strcat('results_number_cog_disconn_',num2str(start_time(1)), '_' , num2str(start_time(2)), '_' , num2str(start_time(3)), '_' , num2str(start_time(4)), '_' , num2str(start_time(5)));
save(savename,'results');

dlmwrite([savename '.edge'], results.resulting_statistics_binary, 'delimiter', '\t', 'precision', 4);

% close the waitbar
close(f)
fprintf('Done')