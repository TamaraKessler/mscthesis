% Pre-Allocate some memory
path  = struct;
hyper = struct;
hypo  = struct;
full_map = struct;

% Set path
path.main = 'D:\Tamara\LesionMaps\to be added\';
path.out  = 'D:\Tamara\LesionMaps\to be added\added';

% RHLM numbers of patients whose maps need to be added
%pnum = [20 28 32 55 119 314 356 476 480 500];
pnum = [493];
% Transform for easier handling
pnum = pnum';

% For every patient whose maps need to be added
for i_pat = 1:length(pnum)
    
    % Set path to patient's subfolder
    path.sub = fullfile(path.main, num2str(pnum(i_pat)));
    folder   = dir([path.sub '\*.nii']);
    
    % Get file names of the maps
    hyper.name = folder(1).name;
    hypo.name  = folder(2).name;
    
    % Get paths to the maps
    hyper.path = fullfile(path.sub, hyper.name);
    hypo.path  = fullfile(path.sub, hypo.name);
    
    % Read in Headers of the maps
    hyper.header = spm_vol(hyper.path);
    hypo.header  = spm_vol(hypo.path);
    
    % Now read in the volumes
    hyper.data = spm_read_vols(hyper.header);
    hypo.data  = spm_read_vols(hypo.header);
    
    % Add the two maps
    full_map.data = hyper.data + hypo.data;
    % If there are any overlaps - aka 2s in the added map
    if ismember(2, full_map.data)
        % Correct for this by replacing the 2s with 1s
        full_map.data(full_map.data==2) = 1;
    end
    
    % Copy over the header from one of the OG maps
    full_map.header = hyper.header;
    % Adjust filename in header 
    new_name = ['RHLM_' num2str(pnum(i_pat)) '_bwshypo_and_hyper_label.nii'];
    full_map.header.fname = fullfile(path.out, new_name);
    full_map.header.private.dat.fname = full_map.header.fname;
    
    % Write image
    spm_write_vol(full_map.header, full_map.data);
    
end
