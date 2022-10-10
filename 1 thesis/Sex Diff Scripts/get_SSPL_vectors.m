% Get idcSSPLs per person

pnum = 103;
mean_SSPL = zeros(pnum,2);
max_SSPL = zeros(pnum,2);

path.female = 'D:\Tamara\Diskonnektionen\SSPLs\female\delta\SSPLs';
folder.female = dir(path.female);
folder.female(1,:) = [];
folder.female(1,:) = [];

path.male = 'D:\Tamara\Diskonnektionen\SSPLs\male\delta\SSPLs';
folder.male = dir(path.male);
folder.male(1,:) = [];
folder.male(1,:) = [];

for i_female = 1:pnum
    
    tmp = load(fullfile(path.female, folder.female(i_female).name)).delta_sspl_matrix;
    mean_SSPL(i_female,1) = mean(tmp(:)); 
    max_SSPL(i_female,1) = max(tmp(:));
    clear tmp
    
end

for i_male = 1:pnum
    
    tmp = load(fullfile(path.male, folder.male(i_male).name)).delta_sspl_matrix;
    mean_SSPL(i_male,2) = mean(tmp(:)); 
    max_SSPL(i_male,2) = max(tmp(:));
    clear tmp
    
end

saveloc = 'D:\Tamara\Diskonnektionen\SSPLs\';
writematrix(mean_SSPL,fullfile(saveloc, 'meanSSPL_vectors_delta.csv'));
writematrix(max_SSPL,fullfile(saveloc, 'maxSSPL_vectors_delta.csv'));
