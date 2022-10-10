%% ************************* MOVE MATFILES *************************** %%
% 
% Written by Tamara Keßler, 06/2022
%
%%

% For convience, move every patient's disconnection map matfiles into the
% same folder

%% Set up

clear; clc;

% Set up paths
path.disconmaps = 'D:\Tamara\LesionMaps\Shifted\all';
path.disconmats = 'D:\Tamara\LesionMaps\Shifted\all\mat_files';
folder = dir([path.disconmaps, '\*.mat']);

n_pat = size(folder,1);

%%

% For every patient
for i_pat = 1:n_pat
    
    % Retrieve the name & path to the patient's mat file
    pat_name = folder(i_pat).name;
    path2mat = fullfile(path.disconmaps, pat_name);
    
    % Move file to according folder
    movefile(path2mat,path.disconmats);
    
    clear pat_name path2mat 
    
end %i_pat