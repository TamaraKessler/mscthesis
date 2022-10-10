%% ************************* ADD LESION MAPS *************************** %%
% 
% Written by Tamara Keßler, 06/2022
%
%%

% For convenience, this script moves all parcel-wise disconnection maps into
% the same folder

%% Set up

clear; clc;

% Set up paths
dir.patnum = 'D:\Tamara\Data\Patientlists_Groups';
dir.disconnection = 'D:\Tamara\Diskonnektionen\LQT_output';
dir.disconmaps = 'D:\Tamara\Diskonnektionen\SSPLs\indirect';

% Read in patient numbers
RHLM_numbers = importdata(fullfile(dir.patnum, 'all_lessmen.txt'));
n_pat = size(RHLM_numbers,1);

%%

fprintf('Moving disconnection files...\n');

% For every patient
for i_pat = 2:n_pat
    
    fprintf('Map %d/%d moved \n',i_pat,n_pat);
    
    % Retrieve name & path of the file
    pat_name = [RHLM_numbers{i_pat} '.nii'];
    path2map = fullfile(dir.disconnection, pat_name, 'Parcel_SSPL');
    map_name = [pat_name, '_BN_Atlas_indirect_SDC.mat'];
    %map_name = [pat_name, '_BN_Atlas_percent_tdi.nii'];
    
    fullname = fullfile(path2map, map_name);
    
    % copy file into new folder
    copyfile(fullname,dir.disconmaps)
    
    clear pat_name path2map map_name full_name
    
end %i_pat

fprintf('Done.\n');