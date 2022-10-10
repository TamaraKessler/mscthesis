clear; clc;

BN = niftiread('C:\Users\Sperber\Desktop\BN\BN_Atlas_246_1mm.nii');

tmp = BN;
tmp(tmp<231) = 0;
thal = double(tmp>0);

thalsize = sum(thal(:));


path.female = 'D:\Tamara\LesionMaps\Shifted\all_female';
folder.female = dir(fullfile(path.female, '*.nii'));
path.male = 'D:\Tamara\LesionMaps\Shifted\all_male';
folder.male = dir(fullfile(path.male, '*.nii'));

pnum = 103;

% read in lesion map
% read in BN

female_thallesion = zeros(pnum,3);
fem_thallesion_pnum = cell(pnum,1);
tempmap = zeros(182,218,182);
maps = zeros(182,218,182,pnum);

for i_female = 1:pnum
    
    file = fullfile(path.female, folder.female(i_female).name);
    name = extractBetween(file,'s_','_b');
    
    tmp = niftiread(file);
    tempmap = BN.*tmp;
    
    tempmap(tempmap < 231) = 0;
    tempmap_bin = double(tempmap>0);
    
    if max(tempmap(:)) == 0
        female_thallesion(i_female,:) = 0;
    else
        female_thallesion(i_female,1) = 1;
        female_thallesion(i_female,2) = sum(tempmap_bin(:));
        female_thallesion(i_female,3) = (sum(tempmap_bin(:))/thalsize)*100;
        fem_thallesion_pnum{i_female,1} = name{1,1};
    end
end

male_thallesion = zeros(pnum,3);
male_thallesion_pnum = cell(pnum,1);
tempmap = zeros(182,218,182);
maps = zeros(182,218,182,pnum);

for i_male = 1:pnum
    
    file = fullfile(path.male, folder.male(i_male).name);
    name = extractBetween(file,'s_','_b');
    
    tmp = niftiread(file);
    tempmap = BN.*tmp;
    
    tempmap(tempmap < 231) = 0;
    tempmap_bin = double(tempmap>0);
    
    if max(tempmap(:)) == 0
        male_thallesion(i_male,:) = 0;
    else
        male_thallesion(i_male,1) = 1;
        male_thallesion(i_male,2) = sum(tempmap_bin(:));
        male_thallesion(i_male,3) = (sum(tempmap_bin(:))/thalsize)*100;
        male_thallesion_pnum{i_male,1} = name{1,1};
    end
end
    
nr_fem_thallesions = sum(female_thallesion(:,1));
mean_fem_thallesion = mean(female_thallesion(:,2));
mean_fem_thaldmg = mean(female_thallesion(:,3));
nr_male_thallesions = sum(male_thallesion(:,1));
mean_male_thallesion = mean(male_thallesion(:,2));
mean_male_thaldmg = mean(male_thallesion(:,3));

fem_thallesion_pnum = fem_thallesion_pnum(~cellfun('isempty',fem_thallesion_pnum));
male_thallesion_pnum = male_thallesion_pnum(~cellfun('isempty',male_thallesion_pnum));

saveloc = 'D:\Tamara\LesionMaps';

writecell(fem_thallesion_pnum,fullfile(saveloc,'pnums_fem_thalamiclesions.csv'));
writecell(male_thallesion_pnum,fullfile(saveloc,'pnums_male_thalamiclesions.csv'));

results = struct;
results.nr_tha_lesions = nr_fem_thallesions + nr_male_thallesions;
results.nr_female_tha_lesions = nr_fem_thallesions;
results.nr_male_tha_lesions = nr_male_thallesions;
results.tha_size_vx = thalsize;
results.mean_female_tha_lesionsize_vx = mean_fem_thallesion;
results.mean_male_tha_lesionsize_vx = mean_male_thallesion;
results.female_overview = female_thallesion;
results.female_pnums_with_tha_lesion = fem_thallesion_pnum;
results.male_overview = male_thallesion;
results.male_pnums_with_tha_lesion = male_thallesion_pnum;

cd(saveloc)
save results

% multiply
% are there voxels >231?
% count