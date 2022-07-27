path.lesions = 'D:\Tamara\LesionMaps\Originals\masked\all_male\mat_files';
path.data = 'D:\Tamara\Data\Patientlists_Groups';
folder = dir([path.lesions, '\*.mat']);

RHLM_numbers = importdata(fullfile(path.data, 'all_male.txt'));

n_pat = size(folder,1);

for i_pat = 1:n_pat
    
    new_name = fullfile(path.lesions, [RHLM_numbers{i_pat} '.mat']);
    old_name = fullfile(path.lesions, folder(i_pat).name);
    
    movefile(old_name,new_name);
end