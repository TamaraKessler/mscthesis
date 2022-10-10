
clear all; clc;

path.wkdir = 'D:\Tamara\LesionMaps\Originals\masked';

group.list = {'all', 'female', 'male'};
group.idx = listdlg('PromptString','Select the corresponding group', ...
    'SelectionMode','single','ListString',group.list);

path.folder = fullfile(path.wkdir, group.list{group.idx}, 'niistat_results');

cd(path.folder);

path.BN = 'C:\Users\Sperber\Desktop\BN';


%%

path.templates = 'C:\Users\Sperber\Desktop\templates';

AAL = niftiread([path.templates '\aal.nii']);
AAL = double(AAL);

tmp = niftiread('threshZresultsmean_z.nii');
tmp = double(tmp>0);

atlasmap.AAL = AAL.*tmp;

ROIs = unique(AAL);

%%

