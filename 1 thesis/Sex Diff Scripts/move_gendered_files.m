%% ************************* MOVE FILES *************************** %%
% 
% Written by Tamara Keßler, 06/2022
%
%%

% For convience, move every patient's disconnection map matfiles into the
% same folder

%%

clear; clc;

% Set up paths
path.in = 'D:\Tamara\Diskonnektionen\SSPLs\all\delta\';
path.data = 'D:\Tamara\Data\Patientlists_Groups\';
folder = dir([path.in, '*.mat']);

path.out = 'D:\Tamara\Diskonnektionen\SSPLs\female\delta';

RHLM_numbers = importdata(fullfile(path.data, 'all_female_nozeros.txt'));

% Allocate some memory
n_pat = size(folder,1);

%%

% For every patient
for i_pat = 1:n_pat
    
    % Retrieve old path & name
    pat_num = extractBefore(folder(i_pat).name, ".nii_BN");
    filename = fullfile(path.in, folder(i_pat).name);
    
    if ismember(pat_num,RHLM_numbers)
    
        % Copy file to final folder
        copyfile(filename, path.out)
    
    end %if
        
end %i_pat

fprintf('Done.\n');