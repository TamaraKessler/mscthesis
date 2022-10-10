%% ********************** MOVE NII ORIGIN FOR LQT ************************ %%
% 
% Written by Tamara Keßler, 06/2022
%
%%

% This script moves our images from the 181x217x181 space used in our group
% to 182x218x182 space for the Lesion Quantification Toolkit

%% Set up

clear; clc;

% Pre-Allocate some memory
orig  = struct;
shift = struct;
path  = struct;
path.main = uigetdir('D:\Tamara\LesionMaps\Originals\masked', 'Choose input folder');
path.out  = uigetdir('D:\Tamara\LesionMaps\Shifted', 'Choose output folder');

fileList = dir(fullfile(path.main, '*.nii'));

% Get number of scans
n_scans = size(fileList,1);

%%

fprintf('++++++++++++++++++++++++\n');

% For every map
for i_map=1:n_scans
    
    fprintf('+++ File %d/%d moved +++\n', i_map, n_scans);
    
    % Get the path to the lesion map
    path.map = fullfile(path.main,fileList(i_map).name);
    % Read in the header
    orig.header = spm_vol(path.map);
    % Read in the map
    orig.img = spm_read_vols(orig.header);
    % Save the filename of the OG image
    orig.name = fileList(i_map).name;

    if ~isequal(orig.header.dim, [181 217 181]) %check if dimensions are equal
        warning('Skipping %s: dimensions not 181x217x181\n', path.scan);
        continue;
    end
    
    % Create a new image, only consisting of zeros
    shift.img = zeros(182,218,182);
    % Paste the actual image onto it
    shift.img(2:182,2:218,2:182) = orig.img; 
    % Flip image
    shift.img = flip(shift.img,1);
    
    % Copy over header from original map
    shift.header = orig.header;
    % Adjust bounding box size etc
    shift.header.mat = [-1 0 0 91;0 1 0 -127;0 0 1 -73; 0 0 0 1]; 
    shift.header.dim = [182 218 182];
    % Change file name
    shift.name = ['s_' orig.name];
    shift.header.fname = fullfile(path.out, shift.name);
    shift.header.private.dat.fname = shift.header.fname;
    
    % Write shifted image
    spm_write_vol(shift.header,shift.img);
    
end % i_map

fprintf('++++++++++++++++++++++++++\n');
fprintf('\n');
fprintf('Finished.\n');