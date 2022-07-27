%% *************************** Mask Ventricles ****************************
%-----------------------------------------------------------------------
% Lisa Röhrig, 01-12-2020
%-----------------------------------------------------------------------
%
% Code template by Chris (2014)
% Applies mask to lesion maps, cuts out ventricles/cerebellum 
% Nifti-mask & nifti-lesions to apply mask on are needed
% (E.g., mask_ch2-template.nii that is based on the ch2-template)

clear; clc;

%path = 'D:\Lisa\Projekt WMH\Lisa WMH_30112020\Data\Modified\';
path.in = uigetdir('D:\Tamara\LesionMaps\Originals\mod', 'Select input folder');
cd(path.in);

path.out = uigetdir('D:\Tamara\LesionMaps\Originals\masked', 'Select output folder');

% Select mask and lesion files
%mask = 'D:\Lisa\Projekt WMH\Lisa WMH_30112020\mask_ch2-template.nii'; % ch2-template
mask = 'C:\Users\Sperber\Desktop\templates\mask_ch2-template.nii';

% if ~exist('mask','var')
%     mask = spm_select(1,'image','Select mask'); 
% end

if ~exist('lesions','var')
    lesions = spm_select(inf,'image','Select lesion-image[s] to apply mask on'); 
end
  
fprintf('Apply Mask on Selected Lesion Maps....\n');

% Number of files
n = length(lesions);

% Load mask file
hdr_mask = spm_vol(mask);    
img_mask = spm_read_vols(hdr_mask);

mask = (img_mask == 0);

% For all lesion files
for i = 1:n
    % Load lesion file
    curr_lesion = deblank(lesions(i,:));
    hdr_lesion = spm_vol(curr_lesion);    
    img_lesion = spm_read_vols(hdr_lesion);

    % Apply mask on lesion map
    img_new = img_lesion .* mask;

    % Write to new matrix (same folder structure required)
    path_new = replace(hdr_lesion.fname, 'Modified', 'Masked');
    filename_new = strcat(extractBetween(path_new, [path.in, '\'], '.nii'), '_masked.nii');
    fullname = fullfile(path.out, filename_new(1));
    hdr_new = hdr_lesion;
    hdr_new.fname = fullname{1};
    spm_write_vol(hdr_new, img_new);
    file = extractAfter(path_new, 'Masked\');
    fprintf('(%d) File \t>>%s<<\t....masked\n', i, file);
end
