%% ******************** DELETE EXCLUDED FILES *********************** %%
% 
% Written by Tamara Keßler, 06/2022
%
%%

% This script deletes files associated with specified patient numbers

%% Set up

clear; clc;

% Set up paths
path.lesions = 'D:\Tamara\LesionMaps\Originals\masked\all\mat_files';
path.data = 'D:\Tamara\Data\Patientlists_Groups';
folder = dir([path.lesions, '\*.mat']);

RHLM_numbers = importdata(fullfile(path.data, 'all.txt'));

% Allocate some memory
n_pat = size(folder,1);

tbd = [452 454 476 478 481 482 485 486 490 492 498 502 503 504 505 510];

%%

% For every patient
for i_pat = 1:n_pat
    
    % Retrieve old path & name
    pat_num = str2double(extractBetween(folder(i_pat).name, "RHLM_", ".mat"));
    filename = fullfile(path.lesions, folder(i_pat).name);
    
    if ismember(pat_num,tbd)
    
        % Rename by "moving" the file from the old name to the new one
        delete(filename)
    
    end %if
        
end %i_pat