% Set paths
path = 'D:\Tamara\LesionMaps\all';
folder = dir([path '\*.nii']);

% Count how many scans there are
n_scans = size(folder,1);

% Pre-Allocate some memory
lesion_ccm = zeros(n_scans,2);

% For every scan
for i_scan = 1:n_scans 

    % Load the nift file
    img_name = folder(i_scan).name;
    fullpath = fullfile(path, img_name);
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
