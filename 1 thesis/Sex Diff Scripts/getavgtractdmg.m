
% Set paths

path.in = 'D:\Tamara\Diskonnektionen\LQT_output';
folder = dir(path.in);

folder(1) = [];
folder(1) = [];

concat_tractmats = zeros(206,70);

% Read in every patient's mat tract disconnection file

for i = 1:length(folder)
    
    tmp.name = folder(i).name;
    path.tractmat = fullfile(path.in, tmp.name, 'Tract_Disconnection', ...
        [tmp.name '_percent_discon_tracts.mat']);
    tmp.strct = open(path.tractmat);
    tmp.tractmat = tmp.strct.tract_discon;

    concat_tractmats(i,:) = tmp.tractmat;
    
    clear tmp
end

% Average

avg_tractmat = mean(concat_tractmats,1);