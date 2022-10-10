%% ************ DISCONNECTION VISUALISATION (GROUP LEVEL) *************** %%
% 
% Written by Tamara Keßler, 06/2022
%
%%

% This script creates group-level whole-brain disconnection visualisations
% (i.e., average disconnection& overlay of disconnections), that can be
% viewed in MRIcron

%% Set up

clear; clc;

% Set up paths
wkdir = "D:\Tamara";
prompt = "Select input folder containing .mat files";
folder = uigetdir(wkdir, prompt); % May only contain the relevant mat files

%% Read data

% Read imaging data
fileList = struct2cell(dir(fullfile(folder, '*.mat')));
temp=load(strcat(fileList{2,1}, '\', fileList{1,1})).rest.dat; %load first file to get matrix size

images_3d = zeros(size(temp,1),size(temp,2),size(temp,3),length(fileList)); %create blank for 2d images

for i=1:length(fileList)
    images_3d(:,:,:,i)=load(strcat(fileList{2,i}, '\', fileList{1,i})).rest.dat;
end

%% create a 2d mask for the imaging data

% then we need to remove rarely/never affected connections
images_3d_binary=double(images_3d>0);% binarise the image, every disconnection >0 is set to 1
image_3d_sum=sum(images_3d_binary,4);

% % you can choose to only test connections that are affected in at least threshold = x
% % patients (I set it to 5, but feel free to change it)
threshold = 5;
image_3d_sum_thresh = image_3d_sum;
image_3d_sum_thresh(image_3d_sum_thresh<threshold)=0;

%%

name = 'male.nii';
saveloc = 'D:\Tamara\Diskonnektionen\Disconnection_Maps\Overlays';
avg = 'Average_Discon_';
overlap = 'Overlap_Discon_';

%% create an average disconnection plot

mean_matrix_3d = mean(images_3d,4);

header = spm_vol('D:\Tamara\Diskonnektionen\Disconnection_Maps\all\RHLM_1.nii_BN_Atlas_percent_tdi.nii');
savename = fullfile(saveloc, [avg, name]);
header.fname = savename;
%header.fname = 'D:\Tamara\DisconnectionMaps\Average_Discon_conM.nii';
header.private.dat.fname = header.fname;

spm_write_vol(header, mean_matrix_3d);

clear savename header

%% create an overlay disconnection plot

matrix_3d_sum = sum(images_3d_binary,4);

header = spm_vol('D:\Tamara\Diskonnektionen\Disconnection_Maps\all\RHLM_1.nii_BN_Atlas_percent_tdi.nii');
savename = fullfile(saveloc, [overlap, name]);
header.fname = savename;
%header.fname = 'D:\Tamara\DisconnectionMaps\Overlay_Discon_conM.nii';
header.private.dat.fname = header.fname;

spm_write_vol(header, matrix_3d_sum);