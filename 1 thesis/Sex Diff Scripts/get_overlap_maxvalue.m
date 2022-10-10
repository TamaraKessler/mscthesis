%path = 'D:\Tamara\LesionMaps\Originals\masked\overlay\subtraction\';
path = 'D:\Tamara\Diskonnektionen\Disconnection_Maps\Overlays\'
folder = dir(path);

n_files = size(folder,1);
maxvals = zeros(n_files,1);

for i_files = 1:n_files
    
    if folder(i_files).bytes == 0
        continue;
    end;
    
    name = folder(i_files).name;
    header = spm_vol([path, name]);
    img = spm_read_vols(header);
    
    maxvals(i_files) = max(img(:));
    
    clear name header img
    
end

maxvals = nonzeros(maxvals);