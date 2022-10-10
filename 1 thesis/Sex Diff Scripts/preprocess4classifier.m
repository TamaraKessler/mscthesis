
%% Set up

clear; clc;

%path.wkdir = 'D:\Tamara\Diskonnektionen\Disconnection_Maps';
path.wkdir = 'D:\Tamara\LesionMaps\Shifted';

prompt.folder = 'Select input folder';
path.folder = uigetdir(path.wkdir, prompt.folder);

cd(path.folder);

group.list = {'all', 'female', 'male'};
group.idx = listdlg('PromptString','Select the corresponding group', ...
    'SelectionMode','single','ListString',group.list);

%% Read in data

fileList = struct2cell(dir(fullfile(path.folder, '*.mat')));
% short check: mat-files need to be in the same order as in the
% behavioural file; ideally, double check this here!
temp=load(strcat(fileList{2,1}, '\', fileList{1,1})).lesion.dat; %load first file to get matrix size
% IMPORTANT: depending on the format of the data you have to modify this
% line. This line works with the LQT output, which is a struct that
% contains the a variable 'pct_sdc_matrix'. The '.' is the syntax to access
% an element within a struct. If you open a text file, an excel file or a
% matlab file without such a struct, you have to adapt this line

matsize = size(temp);

images_3d = zeros(matsize(1),matsize(2),matsize(3),length(fileList)); %create blank for 3d images

for i_file=1:length(fileList)
    images_3d(:,:,:,i_file)=load(strcat(fileList{2,i_file}, '\', fileList{1,i_file})).lesion.dat;
    % IMPORTANT: same thing as in the previous comment!
end

clear i_file

%% Mask

mask=ones(matsize);

images_3d_binary = double(images_3d);
images_3d_binary_sum = sum(images_3d_binary,4);
mask(images_3d_binary_sum < 5) = 0;
mask_vect=mask(:);

%% Vectorise -> 2D

images_vect=zeros(length(fileList),sum(mask(:))); %create blank

for i_file=1:length(fileList)
    temp=images_3d_binary(:,:,:,i_file);
    temp=temp(:);
    images_vect(i_file,:)=temp(mask_vect==1);
end

%% PCA

% want: score + explained (<- cumsum(explained))
[~, score, ~, ~, explained, ~] = pca(images_vect);

% sumvar = 0;
% nr_components = 0;
% 
% while sumvar < 95
%     
%     nr_components = nr_components + 1;
%     sumvar = sumvar+explained(nr_components);
%     
% end

var_thr = 95; % Variance threshold
cum_var = cumsum(explained); % Cumulative sum of variance explained by components
for h = 1:length(cum_var)
    if cum_var(h) > var_thr
        n_comp = h; % Number of components
        fprintf('\nMaintain %d components (total = %d components) with %f cumulative percent explained variance\n\n', n_comp, size(score,2), cum_var(h));
        break
    end
end

score_red = zeros(length(fileList),n_comp);
score_red = score(:,1:n_comp);

%% Plot

%plot(sqrt(sum(score.^2,1)))

% Make plot
t_day = datestr(now, 1);
t_time = strcat(datestr(now, 'HH'), '-', datestr(now, 'MM'));
figure('Name', 'explained Var')
l = size(score,2);
line(1:l, cum_var, 'marker', 's', 'color', 'k', 'markerfacecolor', 'white', 'MarkerSize', 5);
yline(var_thr, '--k'); % Add threshold line
xlabel('Number of principal components');
ylabel('Cumulative variance explained');
legend('Data', 'Threshold', 'Location', 'southeast');
%saveas(gcf, strcat(img_dir, '\Fig_PCA_', title, '_rs_', num2str(r), '_', data_txt, '_', t_day, '_', t_time, '.png'));


%% Scale

minval = min(score_red(:));
maxval = max(score_red(:));
score_red_normed = ((score_red-minval)/(maxval-minval));

tmp = readmatrix('D:\Tamara\LesionMaps\Originals\masked\all\mat_files\renamed\all_behaviour_neg.csv');
behaviour = tmp(:,2);

clear minval maxval
minval = min(behaviour(:));
maxval = max(behaviour(:));

behaviour_normed = ((behaviour-minval)/(maxval-minval));

%%

cd('D:\Tamara\Diskonnektionen\LibSVM\lesionmaps');

clearvars -except score_red_normed behaviour_normed

save data4libsvm
%writematrix(behaviour_normed,'labels.csv');
writematrix(score_red_normed, 'instance_matrix.csv');