function FixOrigin_nii(fnmOK, fnms)
% Function changes origin of nii-image
% -------------------------------------------------------------------------
% Original by Chris Rorden and modified by Lisa Röhrig
% -------------------------------------------------------------------------
%
% **************** => APPLIES FOR IMAGES THAT HAVE ORIGIN [-90, -126, -72] 
% AND SHOULD HAVE ORIGIN [-90, -125, -71] ********************************
% -> changes wrong origin (Bounding Box) in header file + moves image accordingly
%    (images have to be flipped and moved 1 voxel in 2 dimension-directions [y, z])
% -> 1 image with correct origin needed that is used to adapt the other files with a wrong origin
% -> files are skipped if they already have the desired origin
% ####### code has to be modified in case of another desired/available origin

% ---------------------- Function by Chris Rorden ------------------------
% Fix NIfTI images  with incorrect origin/s_form/q_form values
% fnmOK: image with correct origin (other images modified to match this)
% fnms: image(s) to fix
% Note: designed for 3D images, 4D will require modifications
% Examples
% nii_fixorigin(); %use gui
% nii_fixorigin('a.nii','b.nii');
% nii_fixorigin('a.nii',strvcat('b.nii','d.nii'));
 
%clear; clc;
path.in = uigetdir('D:\Tamara\LesionMaps\Originals\OG\', 'Select input folder');
path.out = uigetdir('D:\Tamara\LesionMaps\Originals\mod', 'Select output folder');

if ~exist('fnmOK','var')
    fnmOK = spm_select(1,'image','Select image with correct origin'); 
end
if ~exist('fnms','var')
    fnms = spm_select(inf,'image','Select image[s] to correct'); 
end
hdrOK = spm_vol(fnmOK);
    
for i=1:size(fnms,1)
    fnm = deblank(fnms(i,:));
    hdr = spm_vol(fnm);
    if ~isequal(hdr.dim, hdrOK.dim)
        warning('Skipping %s: dimensions differ\n', fnm);
        continue;
    end
    if ~isequal(hdr.mat, hdrOK.mat) % only if origin is not equal
        img = spm_read_vols(hdr);

    % ------------------------- New part - Flip ---------------------------
        left = img(1:91, :, :); % left side of lesion image matrix
        right = img(92:181, :, :); % right side of lesion image matrix
        ones_left = sum(left, 'all'); 
        ones_right = sum(right,'all');
        if ~isequal(hdr.mat(1, 4), hdrOK.mat(1, 4)) % check if the R/L dimension is to be flipped
           fprintf('\nFlipping image... %s\n', fnm); 
           img = flip(img, 1);
           left_flip = img(1:91,:,:);
           right_flip = img(92:181,:,:);
           ones_left_flip = sum(left_flip, 'all');
           ones_right_flip = sum(right_flip, 'all');
           fprintf('Voxels lesioned left: %s\n', num2str(ones_left_flip)); %to double check
           fprintf('Voxels lesioned right: %s\n', num2str(ones_right_flip));    
        else
           fprintf('\nImage is not to be flipped... %s\n', fnm); 
           fprintf('Voxels lesioned left: %s\n', num2str(ones_left)); %to double check
           fprintf('Voxels lesioned right: %s\n', num2str(ones_right));    
        end
    % -------------------------- New part - Change Image ------------------
        new_img = zeros(181, 217, 181);
        new_img(:, 1:216, 1:180) = img(:, 2:217, 2:181); % move image 1 voxel in y- and z-direction

        % Create new header
        new_hdr = hdr;
        new_hdr.mat = hdrOK.mat;
        new_hdr.private.mat = hdrOK.private.mat;
        new_hdr.private.mat0 = hdrOK.private.mat0;
        filename = strsplit(hdr.fname, '.nii');
        filename = extractAfter(filename, [path.in, '\']);
        new_filename = fullfile(path.out, filename);
        new_hdr.fname = char(strcat(new_filename(1), '_mod.nii'));
        spm_write_vol(new_hdr, new_img);
    end
end
end