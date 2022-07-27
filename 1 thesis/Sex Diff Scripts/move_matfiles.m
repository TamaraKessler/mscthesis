path.patnum = 'D:\Tamara\Data\Patientlists_Groups';
path.disconmaps = 'D:\Tamara\DisconnectionMaps\all';
path.disconmats = 'D:\Tamara\DisconnectionMaps\all\mat_files_FA';
folder = dir([path.disconmaps, '\*.mat']);

RHLM_numbers = importdata(fullfile(path.patnum, 'all.txt'));
n_pat = size(folder,1);

for i_pat = 1:n_pat
    
    pat_name = folder(i_pat).name;
    path2mat = fullfile(path.disconmaps, pat_name);
    %path2map = fullfile(dir.disconnection, pat_name, 'Disconnection_Maps');
    
    %fullname = fullfile(path2map, map_name);
    
    movefile(path2mat,path.disconmats);
    
    clear pat_name path2mat 
    
end %i_pat