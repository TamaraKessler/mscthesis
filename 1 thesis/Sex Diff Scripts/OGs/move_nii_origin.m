% this script moves our images from the 181x217x181 space used in our group
% to 182x218x182 space for the Lesion Quantification Toolkit
clear
clc

% Pre-Allocate some memory
path = struct;
path.main = 'D:\Tamara\LesionMaps\';
path.out  = 'D:\Tamara\ShiftedMaps';

%fileList = struct2cell(dir(fullfile(path.main, '*.nii')));
fileList = dir(fullfile(path.main, '*.nii'));

n_scans = size(fileList,1);

% For every scan
for i_scan=1:n_scans
    
    path.scan = fullfile(path.main,fileList(i_scan).name);
    header=spm_vol(path{1,1});
    img=spm_read_vols(header);

    if ~isequal(header.dim, [181 217 181]) %check if dimensions are equal
        warning('Skipping %s: dimensions not 181x217x181\n', path{1});
        continue;
    end
    
    new_img=zeros(182,218,182);
    new_img(2:182,2:218,2:182)=img; 
    % Flip image
    new_img = flip(new_img,1);
    
    new_header=header;
    new_header.mat=[-1 0 0 91;0 1 0 -127;0 0 1 -73; 0 0 0 1]; 
    new_header.dim=[182 218 182];
    new_header.fname=char(strcat(['D:\Lisa\Projekt WMH-Diskonnektionen\Data\LesionMaps\' save_txt],'t_', fileList(1,i)));
    spm_write_vol(new_header,new_img);
    
end % i_scan