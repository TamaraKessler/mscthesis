% Pre-Allocate some memory
path = struct;
base = struct;
mid  = struct;
top  = struct;
full = struct;

% Set paths
path.main = 'D:\Tamara\Scans\add';
path.out  = 'D:\Tamara\Scans\add\added';

% RHLM numbers of patients whose scans need to be concatenated
pnum = [292 295 296 323 333];
% Transform for easier handling
pnum = pnum';

% For every patient whose maps need to be added
for i_pat = 1:length(pnum)
    
    % Set path to patient's subfolder
    path.sub = fullfile(path.main, num2str(pnum(i_pat)));
    folder   = dir([path.sub '\*.nii']);
    
    % Get file names of the scans
    base.name = folder(1).name;
    top.name  = folder(end).name;
    
    % Get paths to the scans
    base.path = fullfile(path.sub, base.name);
    top.path  = fullfile(path.sub, top.name);
    
    % Read in Headers of the scans
    base.header = spm_vol(base.path);
    top.header  = spm_vol(top.path);
    
    % Now read in the volumes
    base.data = spm_read_vols(base.header);
    top.data  = spm_read_vols(top.header);
    
    % Get number of slices
    base.dim = base.header.dim(:,3);
    top.dim  = top.header.dim(:,3);
    
    % First patient's scan is split into thirds
    if i_pat == 1
        mid.name   = folder(2).name;
        mid.path   = fullfile(path.sub, mid.name);
        mid.header = spm_vol(mid.path);
        mid.data   = spm_read_vols(mid.header);
        mid.dim    = mid.header.dim(:,3);
        % Calculate total number of slices
        n_slices   = base.dim + mid.dim + top.dim;
        % Concatenate scans
        full.data  = cat(3,base.data,mid.data,top.data);
    else 
        n_slices   = base.dim + top.dim;
        full.data  = cat(3,base.data,top.data);
    end
    
    newname = ['RHLM_' num2str(pnum(i_pat)) '_full.nii'];
    
    % Create the header
    full.header = top.header;
    full.header.fname = fullfile(path.out, newname);
    full.header.private.dat.fname = full.header.fname;
    full.header.dim = [512 512 n_slices];
    
    % Write image
    spm_write_vol(full.header, full.data);
  
end % i_pat
