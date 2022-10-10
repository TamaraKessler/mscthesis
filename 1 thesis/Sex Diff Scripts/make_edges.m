function make_edges (n,bfrule)

%path = 'D:\Tamara\Diskonnektionen\SSPLs\all\indirect\BF_results\45_BF_nolog.csv';
path = 'D:\Tamara\Diskonnektionen\Parcel_Disconnections\female\BF_results\20_BF_nolog.csv';

tmp = readtable(path);
tmp(:,1) = [];
results_bf = table2array(tmp); clear tmp

%results_bf = readmatrix(['//neurologie/homes/Research_groups/ssmaczny/Desktop/Smaczny 2021_left_ang_disc/Smaczny et al 2021, NeuroImage/Mendeley Data Facts and Figures_updated/results_2d_N',num2str(n),'_BF_',bfrule,'.csv']);
%results_bf(:,1) = [];


tmp = extractBetween(path, 'BF_results\', '_BF');
n = tmp{1}; clear tmp

tmp = extractBetween(path, [n '_BF_'], '.csv');
bfrule = tmp{1}; clear tmp

dlmwrite(['D:\Tamara\Diskonnektionen\Parcel_Disconnections\female\BF_results\',num2str(n),'_BF_',bfrule,'.edge'], results_bf, 'delimiter', '\t', 'precision', 4);

end