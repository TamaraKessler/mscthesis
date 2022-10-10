%% ******************** RENAME FILES FOR NIISTAT *********************** %%
% 
% Written by Tamara Keßler, 07/2022
%
%%

% The NiiStat toolbox requires all files to have the same name as are
% listed in the behavioural data file. This script renames all files
% according to their RHLM patient number

%% Set up

clear; clc;

% Set up paths
path.lesions = 'D:\Tamara\Diskonnektionen\Disconnection_Maps\all\mat_files\renamed';
path.data = 'D:\Tamara\Data\Patientlists_Groups';
folder = dir([path.lesions, '\*.mat']);

%RHLM_numbers = importdata(fullfile(path.data, 'all_lessmen.txt'));

% Allocate some memory
n_pat = size(folder,1);

%%

% For every patient
for i_pat = 1:n_pat
    
    % Retrieve old name
    old_name = fullfile(path.lesions, folder(i_pat).name);
    new_name = [extractBefore(old_name, '.nii_BN') '.mat'];
    
    % Rename by "moving" the file from the old name to the new one
    movefile(old_name,new_name);
    
end %i_pat