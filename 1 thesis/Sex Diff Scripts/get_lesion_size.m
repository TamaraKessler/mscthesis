%% ************************* GET LESION SIZE *************************** %%
% 
% Written by Tamara Keßler, 05/2022
%
%%

% This scripts calculates the lesion size in ccm based on the provided
% binary lesion map (NII) and saves an xlsx file of the results

%% Set up

clear; clc;

% Set paths
path.infolder = 'D:\Tamara\LesionMaps\Originals\OG\all';
folder = dir([path.infolder '\*.nii']);

path.output = 'D:\Tamara\Data';

% Count how many scans there are
n_scans = size(folder,1);

% Pre-Allocate some memory
lesion_ccm = zeros(n_scans,2);

%%

fprintf('Calculating...\n');

% For every scan
for i_scan = 1:n_scans 

    fprintf('Currently at patient %d/%d\n',i_pat,n_pat);
    
    % Load the nifti file
    img_name = folder(i_scan).name;
    fullpath = fullfile(path.infolder, img_name);
    img = niftiread(fullpath);
    % Sum over lesion map to get lesion size in mm
    cmm = sum(img, 'all');
    % Convert cmm to ccm
    ccm = cmm/1000;
    
    % Save lesion size in matrix
    lesion_ccm(i_scan,2) = ccm;
    
    % Extract RHLM number for convience
    RHLM = extractBetween(folder(i_scan).name, 'RHLM_', '_b');
    tmp = convertCharsToStrings(RHLM{1});
    lesion_ccm(i_scan,1) = tmp;
    
    % Clean up workspace
    clear img_name fullpath img cmm ccm RHLM tmp

end % i_scan

cd(path.output);
xlswrite('lesionsize_all.xlsx', lesion_ccm);

fprintf('Done.\n');