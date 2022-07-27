dir.patnum = 'D:\Tamara\Data\Patientlists_Groups';
dir.disconnection = 'D:\Tamara\Diskonnektionen';

RHLM_numbers = importdata(fullfile(dir.patnum, 'all.txt'));
n_pat = size(RHLM_numbers,1);

for i_pat = 1:n_pat
    
    pat_name = RHLM_numbers{i_pat};
    path2map = fullfile(dir.disconnection, pat_name, 'Disconnection_Maps');
    map_name = [pat_name, '_Yeo_100_7_percent_tdi.nii'];
    
    fullname = fullfile(path2map, map_name);
    
    copyfile(fullname,dir.disconmaps)
    
    clear pat_name path2map map_name full_name
    
end %i_pat