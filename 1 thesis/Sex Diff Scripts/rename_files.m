path.lesions = 'D:\Tamara\LesionMaps\Shifted\all';
path.data = 'D:\Tamara\Data\Patientlists_Groups';
folder = dir([path.lesions, '\*.nii']);

RHLM_numbers = importdata(fullfile(path.data, 'all.txt'));

n_pat = size(folder,1);

for i_pat = 1:n_pat
    
    new_name = fullfile(path.lesions, [RHLM_numbers{i_pat} '.nii']);
    old_name = fullfile(path.lesions, folder(i_pat).name);
    
    movefile(old_name,new_name);
end