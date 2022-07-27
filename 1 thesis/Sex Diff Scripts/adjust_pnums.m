path.files = 'D:\Tamara\DisconnectionMaps\Parcel_Disconnection';
%path.data = 'D:\Tamara\Data\Patientlists_Groups';
folder = dir([path.files, '\*.mat']);
path.out = 'D:\Tamara\DisconnectionMaps\Parcel_Disconnection\renamed';

%RHLM_numbers = importdata(fullfile(path.data, 'all.txt'));

n_pat = size(folder,1);

for i_pat = 1:n_pat
    
    old_name = fullfile(path.files, folder(i_pat).name);
    dest = fullfile(path.out, folder(i_pat).name);
    
    copyfile(old_name, dest);
    
    pnum_old = extractBetween(old_name, 'RHLM_', '_Yeo');
    pnum_new = str2num(pnum_old{1})+1000;
    
    new_name = fullfile(path.out, ['RHLM_', num2str(pnum_new), '_parcel.mat']);
    
    movefile(dest,new_name);
end