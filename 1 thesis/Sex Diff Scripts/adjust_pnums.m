%% ********************** ADJUST PATIENT NUMBERS ************************ %%
% 
% Written by Tamara Keßler, 06/2022
%
%%

% The ROI-to-ROI GLM script needs the files to be in the same order as the 
% listed in the behavioural data xlsx file - but Matlab is weird when it
% comes to listing files. To avoid problems, use this file to adjust the
% patient numbers, so they're in the same order

%% Set up

clear; clc;

% Set up paths
path.files = 'D:\Tamara\Diskonnektionen\Parcel_Disconnections';
folder = dir([path.files, '\*.mat']);
path.out = 'D:\Tamara\Diskonnektionen\Parcel_Disconnections\renamed';

n_pat = size(folder,1);

%%

% For every patient
for i_pat = 1:n_pat
    
    % Retrieve filename
    old_name = fullfile(path.files, folder(i_pat).name);
    % Specify output folder
    dest = fullfile(path.out, folder(i_pat).name);
    
    % Copy file into the folder
    copyfile(old_name, dest);
    
    % Extract RHLM number from filename
    pnum_old = extractBetween(old_name, 'RHLM_', '.nii');
    % Add 1000
    pnum_new = str2num(pnum_old{1})+1000;
    
    % Create new name
    new_name = fullfile(path.out, ['RHLM_', num2str(pnum_new), '_parcel.mat']);
    
    % Rename the file my moving it to the new name
    movefile(dest,new_name);
end