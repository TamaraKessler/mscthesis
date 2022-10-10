
%% Set up via User Input

clear; clc;

path.wkdir = 'D:\Tamara\Diskonnektionen\SSPLs';

prompt.folder = 'Select input folder';
path.folder = uigetdir(path.wkdir, prompt.folder);

cd(path.folder);

group.list = {'all', 'female', 'male'};
group.idx = listdlg('PromptString','Select the corresponding group', ...
    'SelectionMode','single','ListString',group.list);

type.list = {'del', 'idc'};
type.idx = listdlg('PromptString','Select the desired SSPL type', ...
    'SelectionMode','single','ListString',type.list);
if type.list{type.idx} == 'del'
    suffix = 'deltaSSPL';
elseif type.list{type.idx} == 'idc'
    suffix = 'idcSSPL';
else
    error('Invalid SSPL type');
end

prompt.behav = 'Select .xlsx file containing behavioural data';
[behav.name, behav.path] = uigetfile('.xlsx', prompt.behav);
behaviour =  xlsread(fullfile(behav.path, behav.name)); clear behav

%% Read in data

fprintf('Reading data...\n');

fileList = struct2cell(dir(fullfile(path.folder, '*.mat')));

if type.list{type.idx} == 'idc'
    
    temp=load(strcat(fileList{2,1}, '\', fileList{1,1})).idc_matrix; %load first file to get matrix size
    images_2d = zeros(length(temp),length(temp),length(fileList)); %create blank for 2d images
    
    for i=1:length(fileList)
        images_2d(:,:,i)=load(strcat(fileList{2,i}, '\', fileList{1,i})).idc_matrix;
    end
    
else 
    
    temp=load(strcat(fileList{2,1}, '\', fileList{1,1})).delta_sspl_matrix; %load first file to get matrix size
    images_2d = zeros(length(temp),length(temp),length(fileList)); %create blank for 2d images
    
    for i=1:length(fileList)
        images_2d(:,:,i)=load(strcat(fileList{2,i}, '\', fileList{1,i})).delta_sspl_matrix;
      
    end
    
end

clear i prompt

%% create a 2d mask for the imaging data

fprintf('Masking data...\n');

mat_size = length(temp);

mask=zeros(mat_size);
% the matrix is a symmetric matrix, so first remove all elements in the
% diagonal and below
for i=1:mat_size
    for j=1:mat_size
        if i<j
            mask(i,j)=1;
        end
    end
end

clear i j

% then we need to remove rarely/never affected connections
images_2d_binary=double(images_2d>0);% binarise the image, every disconnection >0 is set to 1
image_2d_sum=sum(images_2d_binary,3);
image_2d_sum=image_2d_sum.*mask;

threshold = 5;
mask(image_2d_sum<threshold)=0;
mask_vect=mask(:);

% vectorise the 2d images
images_vect=zeros(length(fileList),sum(mask(:))); %create blank
for i=1:length(fileList)
    temp=images_2d(:,:,i);
    temp=temp(:);
    images_vect(i,:)=temp(mask_vect==1);
end

clear i

%% do correlation

fprintf('Doing original correlation...\n');

% allocate some memory
stat_vect_orig=zeros(sum(mask(:)),1);
cors_orig = zeros(length(stat_vect_orig),1);
cors_p_orig = zeros(length(cors_orig),1);

% do correlation with original data

for i = 1:length(stat_vect_orig)
    
    [R,P] = corr(behaviour, images_vect(:,i),'Type','Spearman');
    cors_orig(i) = R;
    cors_p_orig(i) = P;
    
end

clear R P

%% permutation correlation

fprintf('Starting permutations...\n');

defaultStream = RandStream('mlfg6331_64');

perms = 500;
permutation_max_cors=zeros(perms,1); 
cors_vect_perm_temp=zeros(sum(mask(:)),1);

for i_perm = 1:perms
    
    defaultStream.Substream = i_perm;
    rand_behaviour=behaviour(randperm(defaultStream,length(behaviour)));
    
    for i=1:length(stat_vect_orig)
        
        R = corr(rand_behaviour, images_vect(:,i),'Type','Spearman');
        cors_vect_perm_temp(i) = R;
%         cors_p_vect_perm_temp(i) = P(2);
        
    end
    
    tmp = min(cors_vect_perm_temp);
    permutation_max_cors(i_perm) = tmp;
    
end

permutation_max_cors = sort(permutation_max_cors, 'descend');

clear R

%% get results & save them

fprintf('Saving results...\n');

results_2d_cors = zeros(length(mask));

temp=mask(:);
temp(temp==1)=cors_orig;
results_2d_cors = reshape(temp,[length(mask),length(mask)]);

p_thresh = [0.95, 0.99, 0.995, 0.999];

for i_thresh = 1:length(p_thresh)

    %results_2d_thresholded = double(results_2d_cors_p>p_thresh(i_thresh));
    
    permutation_threshold = permutation_max_cors(round(perms*p_thresh(i_thresh)));
    results_2d_thresholded = double(results_2d_cors<permutation_threshold);
    
    results = struct;
    results.p_level = 1-p_thresh(i_thresh);
    results.perm_thresh = permutation_threshold;
    results.input_data = fileList;
    results.time_finished = clock;
    results.mask_tested_idc_SSPLs = mask;
    results.cors_2d = results_2d_cors;
    %results.cors_p_vals = results_2d_cors_p;
    results.sig_cors_thresholded = results_2d_thresholded;
    
    saveloc = fullfile(path.folder, 'corr_results');
    savename = [group.list{group.idx} '_correlation_' suffix '_behaviour_p_' num2str(1000-(p_thresh(i_thresh)*1000))];
    fullname = fullfile(saveloc, savename);
    
    save(fullname,'results');
    
    dlmwrite([fullname '.edge'], results.sig_cors_thresholded, 'delimiter', '\t', 'precision', 4);
    
end

fprintf('Done.\n');

%%
