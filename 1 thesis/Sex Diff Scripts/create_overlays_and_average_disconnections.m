%% read imaging data

wkdir = "D:\Tamara";
prompt = "Select input folder";
folder = uigetdir(wkdir, prompt);

%folder = 'D:\Tamara\DisconnectionMaps\Disconnection_Maps\all\mat_files_rest\renamed'; % this folder contains all relevant disconnection .mat files and NO other .mat files
fileList = struct2cell(dir(fullfile(folder, '*.mat')));
temp=load(strcat(fileList{2,1}, '\', fileList{1,1})).rest.dat; %load first file to get matrix size

images_3d = zeros(size(temp,1),size(temp,2),size(temp,3),length(fileList)); %create blank for 2d images

for i=1:length(fileList)
    images_3d(:,:,:,i)=load(strcat(fileList{2,i}, '\', fileList{1,i})).rest.dat;
end

%% create a 2d mask for the imaging data
% mask=zeros(length(temp));
% % the matrix is a symmetric matrix, so first remove all elements in the
% % diagonal and below
% for i=1:length(temp)
%     for j=1:length(temp)
%         if i<j
%             mask(i,j)=1;images_3d_binary
%         end
%     end
% end

% then we need to remove rarely/never affected connections
images_3d_binary=double(images_3d>0);% binarise the image, every disconnection >0 is set to 1
image_3d_sum=sum(images_3d_binary,4);
% image_2d_sum=image_2d_sum.*mask;

% % you can choose to only test connections that are affected in at least threshold = x
% % patients (I set it to 5, but feel free to change it)
threshold = 5;
image_3d_sum_thresh = image_3d_sum;
image_3d_sum_thresh(image_3d_sum_thresh<threshold)=0;
% mask_vect=mask(:);

% get an overview on the number of patients with a disconnection for each
% voxel
% temp=sort(image_3d_sum_thresh(:),'descend');
% % in 'temp', you can check how many connections are affected in at least x patients
% 
% temp_sparse = sparse(temp);

%% create an average disconnection plot

mean_matrix_3d = mean(images_3d,4);

header = spm_vol('D:\Tamara\DisconnectionMaps\Disconnection_Maps\all\RHLM_1_disconmap.nii');
header.fname = 'D:\Tamara\DisconnectionMaps\Average_Discon_allM.nii';
header.private.dat.fname = header.fname;

spm_write_vol(header, mean_matrix_3d);

%% create an overlay disconnection plot

matrix_3d_sum = sum(images_3d_binary,4);

header = spm_vol('D:\Tamara\DisconnectionMaps\Disconnection_Maps\all\RHLM_1_disconmap.nii');
header.fname = 'D:\Tamara\DisconnectionMaps\Overlay_Discon_allM.nii';
header.private.dat.fname = header.fname;

spm_write_vol(header, matrix_3d_sum);