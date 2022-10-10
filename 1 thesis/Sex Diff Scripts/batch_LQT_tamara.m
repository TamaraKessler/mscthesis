%% LQT batch script
% Lisa Röhrig, modified by Tamara Keßler, 01-2022
clc
clear all

% Add paths to support tools and core functions
addpath('C:\Users\Chris\Lesion_Quantification_Toolkit\Functions');
addpath(genpath('C:\Users\Chris\Lesion_Quantification_Toolkit\Support_Tools'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set up cfg structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Assign relevant paths %%%%%%
% Path to DSI_Studio program
cfg.dsi_path = 'C:\Users\Chris\Desktop\dsi_studio_64\dsi_studio.exe';
% Path to HCP842 .fib template and tractography atlas
cfg.source_path = 'C:\Users\Chris\Lesion_Quantification_Toolkit\Support_Tools\Tractography_Atlas';
% Path to output directory (patient and atlas result directories will be created within the output directory)
cfg.out_path = 'C:\Users\Chris\Desktop\Test_Tamara\output_LQT';
% Path to lesion (pre-registered to MNI template)
cfg.lesion_path = 'C:\Users\Chris\Desktop\Test_Tamara\lesion_maps';
% Path to parcellation (should have identical dimensions to lesion and be in MNI template space)
cfg.parcel_path = 'C:\Users\Chris\Desktop\Test_Tamara\BN_Atlas_246_1mm.nii';

%%%%% Output Filename Options %%%%%%
% Patient ID (used as prefix for output files)
cfg.pat_id = []; % could be a list, but in this example it is being taken from the file names selected later
% File suffix -- used as suffix for output files. Atlas name is recommended (e.g. AAL, Power, Gordon, etc.).
cfg.file_suffix = 'BN_Atlas';

%%%%% Connectivity File Output Options %%%%%%
load('C:\Users\Chris\Lesion_Quantification_Toolkit\Support_Tools\Parcellations\Schaefer_Yeo\Plus_Subcort\Schaefer2018_100Parcels_7Networks_order_plus_subcort.mat');
%cfg.node_label = t.RegionName; %t.RegionName; % n_regions by 1 cell array of strings corresponding to node labels (i.e. parcel names)
%cfg.node_color = t.NetworkID; %t.NetworkID; % n_regions-by-1 array of integer values corresponding to e.g. network assignments or partitions (used to color nodes in external viewers)
%cfg.parcel_coords = [t.X, t.Y, t.Z]; % Parcel coordinates. Used for plotting ball and stick brain graph. If not supplied, they will be estimated from the parcel file

%%%%% Connectivity Options %%%%%%
% connectivity type ('end' or 'pass'): if 'end', connections are defined based on streamline endpoints. If 'pass', connections are defined based on streamline pass-throughs. 'End' is recommended (most conservative).
cfg.con_type = 'end';
% Percent spared threshold for computing SSPLs (e.g. 100 means that only fully spared regions will be included in SSPL calculation; 1 means that all regions with at least 1% spared will be included. Default is 50)
cfg.sspl_spared_thresh = 50;
% smoothing for voxel maps (FWHM in voxel units)
cfg.smooth = 2; 

%%%%%% Navigate to directory containing lesion files %%%%%%
lesion_dir = cfg.lesion_path;
cd(lesion_dir);
lesion_files = dir('*.nii');
% Loop through lesion files and create measures
for i = 1:length(lesion_files)
    % Get patient lesion file and patient ID
    cfg.lesion_path = fullfile(lesion_dir, lesion_files(i).name); % set lesion path to the current lesion
    cfg.pat_id = lesion_files(i).name(1:end); % extract just the patient ID portion of the filename
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Create Damage and Disconnection Measures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get parcel damage for patient
    util_get_parcel_damage(cfg);
    % Get tract SDC for patient
    util_get_tract_discon(cfg);
    % Get parcel SDC and SSPL measures for patient
    util_get_parcel_cons(cfg);
end
% navigate to output directory
cd(cfg.out_path);