%% ************ EVALUATE PARCEL-WISE DISCONNECTIONS ************** %%
% 
% Written by Tamara Keßler, 07/2022
%
%%

%% Set up

clear; clc;

% Set working directory
wkdir = "D:\Tamara\Diskonnektionen";
cd(wkdir);

% Load GLM results
prompt.infolder = "Select input folder";
path.infolder = uigetdir(wkdir, prompt.infolder);
infolder = dir([path.infolder, '\*.mat']);
n_files = size(infolder,1);

cd(path.infolder);

if ~exist('eval_parcel', 'dir')
    mkdir('eval_parcel');
    cd([path.infolder, '\eval_parcel']);
    mkdir('FWER');
    mkdir('FDR');
end;

% [filename.glm, path.glm] = uigetfile('*.mat');
% temp = load(fullfile(path.glm, filename.glm));
% GLM_results = temp.results; clear temp

% User input
group.list = {'all', 'female', 'male'};
group.idx = listdlg('PromptString','Select the corresponding group', ...
    'SelectionMode','single','ListString',group.list);

% Load BN atlas
path.atlas = "C:\Users\Sperber\Desktop\BN";
filename.atlas = "BNA_subregions_edited.xlsx";
atlas = readtable(fullfile(path.atlas, filename.atlas));

corr.list = {'FWER', 'FDR'};
corr.idx = listdlg('PromptString','Select the desired correction type', ...
    'SelectionMode','single','ListString',corr.list);

clear prompt

for i_file = 1:n_files
    
    filename.disconn = infolder(i_file).name;
    
    temp = load(fullfile(path.infolder, filename.disconn));
    GLM_results = temp.results;
    clear temp

%%
    
%     % FWER correction (standard)
%     sparse_disconns = sparse(GLM_results.resulting_statistics_binary);
%     nr_sig_disconns = sum(GLM_results.resulting_statistics_binary(:));

    if corr.list{corr.idx} == "FWER"
        sparse_disconns = sparse(GLM_results.resulting_statistics_binary);
        nr_sig_disconns = sum(GLM_results.resulting_statistics_binary(:));
    elseif corr.list{corr.idx} == "FDR"
        sparse_disconns = sparse(GLM_results.fdr_controlled_binary);
        nr_sig_disconns = sum(GLM_results.fdr_controlled_binary(:));
    else
        error("Invalid correction type");
    end

    if nr_sig_disconns == 0
        continue
    end
    
    spy(sparse_disconns)
    text.title = append("Region-to-Region Disconnections in ", group.list{group.idx}, " patients");
    text.subtitle = append(num2str(nr_sig_disconns), " significant disconnections (p<", num2str(GLM_results.p_level),")");
    title({text.title; text.subtitle});
    xlabel("Parcel Nr.");
    ylabel("Parcel Nr.");

    %%

    % Find all significant (= non-zero) disconnections and get their
    % "coordinates" within the matrix
    [row, col] = find(sparse_disconns);

    %%

    % Count how many significant disconnections every parcel has
    parcels2count = [row; col];
    [counts, parcel_num] = groupcounts(parcels2count);

    %
    n_parcel = length(parcel_num);
    disconn_hubs_ext = cell(n_parcel,6);
    disconn_hubs_min = zeros(n_parcel,2);

    % For every parcel with significant disconnections, count number of
    % disconnections and retrieve anatomical information
    for i_parcel = 1:n_parcel

        % Get index for atlas; needed because always two labels per row
        idx = round(parcel_num(i_parcel)/2);

        % Save relevant info to cell array, i.e., parcel number, anatomical
        % label and number of region-to-region disconnections incl. this parcel
        disconn_hubs_ext{i_parcel,1} = parcel_num(i_parcel);
        disconn_hubs_ext{i_parcel,2} = table2array(atlas(idx,'Lobe'));
        disconn_hubs_ext{i_parcel,3} = table2array(atlas(idx,'Gyrus'));
        disconn_hubs_ext{i_parcel,4} = table2array(atlas(idx,'LeftAndRightHemisphere'));
        disconn_hubs_ext{i_parcel,5} = table2array(atlas(idx,'Var6'));
        disconn_hubs_ext{i_parcel,6} = counts(i_parcel);

        % Save reduced form, i.e., only parcel number & nr. of disconnections
        disconn_hubs_min(i_parcel,1) = parcel_num(i_parcel);
        disconn_hubs_min(i_parcel,2) = counts(i_parcel);


    end %i_parcel

    clear idx i_parcel parcels2count

    %%

    [nr,idx] = max(disconn_hubs_min(:,2));
    ana = append('Parcel ', num2str(disconn_hubs_ext{idx,1}), ' / ', ...
        string(disconn_hubs_ext{idx,3}),'(',string(disconn_hubs_ext{idx,3}),')');
    max_disconn_node = append(ana, ' -> ', num2str(nr), ' Disconnections');

    clear nr idx ana

    %%


    dc_counts = disconn_hubs_ext(:,6);

    disconn_lobes = string(disconn_hubs_ext(:,2));
    cat_lobes  = categorical(string(disconn_lobes));
    lobe_names = categories(cat_lobes);
    lobe_counts = countcats(cat_lobes);

    lobe_disconn_count = zeros(length(lobe_names),1);

    for i_dc = 1:length(disconn_lobes)
        for j_cat = 1:length(lobe_names)
            
            if disconn_lobes(i_dc) == lobe_names(j_cat)

                lobe_disconn_count(j_cat) = lobe_disconn_count(j_cat)+dc_counts{i_dc};

            end
        end
    end

    disconn_hubs_lobes = cell(length(lobe_names),2);
    disconn_hubs_lobes(:,1) = lobe_names;
    disconn_hubs_lobes(:,2) = num2cell(lobe_disconn_count);

    %%

    roi_disconns = cell(length(col),5);

    for i_disconn = 1:length(col)

        idx1 = round(row(i_disconn)/2);
        idx2 = round(col(i_disconn)/2);

        roi_disconns{i_disconn,1} = append(num2str(row(i_disconn)),'-',num2str(col(i_disconn)));

        temp1 = table2array(atlas(idx1,'LeftAndRightHemisphere'));
        temp2 = table2array(atlas(idx1,'Var6'));
        roi_disconns{i_disconn,2} = append(temp1,' / ',temp2);
        clear temp1 temp2;

        temp1 = table2array(atlas(idx2,'LeftAndRightHemisphere'));
        temp2 = table2array(atlas(idx2,'Var6'));
        roi_disconns{i_disconn,3} = append(temp1,' / ',temp2);
        clear temp1 temp2;

        hem1 = mod(row(i_disconn),2);
        hem2 = mod(col(i_disconn),2);

        if hem1~=hem2
            roi_disconns{i_disconn,4} = 'inter';
            roi_disconns{i_disconn,5} = '---';
        else
            roi_disconns{i_disconn,4} = 'intra';
            
            if mod(hem1,2)==0
                roi_disconns{i_disconn,5} = 'right';
            else
                roi_disconns{i_disconn,5} = 'left';
            end
                
        end

    end % i_disconn

    hems      = roi_disconns(:,4);
    cat_hems  = categorical(hems);
    cat_names = categories(cat_hems);
    [mx, idx_max] = max(countcats(cat_hems));
    [mn, idx_min] = min(countcats(cat_hems));

    max_hem_disconns = append(cat_names{idx_max}, 'hemispheric - ', num2str(mx), ' disconnections');
    min_hem_disconns = append(cat_names{idx_min}, 'hemispheric - ', num2str(mn), ' disconnections');

    LR_intra = roi_disconns(:,5);
    cat_LR = categorical(LR_intra);
    cat_LR_names = categories(cat_LR);
    LR_counts = countcats(cat_LR);
    
    counts_left  = LR_counts(cat_LR_names=="left");
    counts_right = LR_counts(cat_LR_names=="right");
    counts_intra = counts_left + counts_right; clear counts_right
    
    LR_percent = counts_left/counts_intra;
    
    left_intra_dcs = append(num2str(round(LR_percent*100)), '% left intra DCs (', ...
        num2str(counts_left), '/', num2str(counts_intra), ')');
    
    clear i_disconn idx1 idx2 hem1 hem2 hems cat_hems cat_names mx mn ...
        idx_max idx_min LR_intra cat_LR cat_LR_names LR_counts counts_left ...
        counts_intra LR_percent
    
    %%
    
    max_nr = 5;
    
    tmp = GLM_results.resulting_glm_t_statistic;
    [maxvals, maxidxs] = maxk(tmp(:),max_nr);
    
    max_disconns = cell(max_nr,4);
    max_row = zeros(max_nr,1);
    max_col = zeros(max_nr,1);
    
    
    for i_disconn = 1:max_nr
        
        % Get matrix indices of max vals
        [max_row(i_disconn), max_col(i_disconn)] = find(tmp == maxvals(i_disconn));
        
        idx1 = round(max_row(i_disconn)/2);
        idx2 = round(max_col(i_disconn)/2);
        
        
        tmp1 = table2array(atlas(idx1, 'Gyrus'));
        tmp2 = table2array(atlas(idx1, 'Var6'));
        max_disconns{i_disconn,1} = append(tmp1,' / ', tmp2);
        clear tmp1 tmp2
        
        tmp1 = table2array(atlas(idx2, 'Gyrus'));
        tmp2 = table2array(atlas(idx2, 'Var6'));
        max_disconns{i_disconn,2} = append(tmp1,' / ', tmp2);
        clear tmp1 tmp2
        
        hem1 = mod(max_row(i_disconn),2);
        hem2 = mod(max_col(i_disconn),2);

        if hem1~=hem2
            max_disconns{i_disconn,3} = 'inter';
        else
            max_disconns{i_disconn,3} = 'intra';
        end
        
        max_disconns{i_disconn,4} = maxvals(i_disconn);
        
    end
    
  

    %%

    results = struct;
    results.perm_num = GLM_results.perm_num;
    results.p_value = GLM_results.p_level;
    results.used_corr = corr.list{corr.idx};
    results.nr_sig_disconns = nr_sig_disconns;
    results.sparse_sig_disconns = sparse_disconns;
    results.roi_disconns = roi_disconns;
    results.max_hemi_dmg = max_hem_disconns;
    results.min_hemi_dmg = min_hem_disconns;
    results.intra_dcs = left_intra_dcs;
    results.disconn_hubs = disconn_hubs_ext;
    results.max_disconn_hub = max_disconn_node;
    results.disconn_lobes = disconn_hubs_lobes;
    results.max_stat_disconns = max_disconns;
    results.time_finished = clock;

    %%

    % User input: Select where results should be saved
    
    path.outfolder = fullfile(path.infolder, 'eval_parcel', corr.list{corr.idx});
    
%     prompt.outfolder = "Select output folder";
%     path.outfolder = uigetdir(wkdir, prompt.outfolder);

    % % User input: Provide results file name
    % prompt.name = "Please provide a name for the results file.";
    % prompt.title = "How do you want to name the file?";
    % dims = [1 60];
    % default = "results";
    %
    % answer = inputdlg(prompt.name, prompt.title, dims, default);
    % filename.output = answer{1};

    p_value  = num2str(10000-(results.p_value*10000));

    savetime = date;
    savename = append('results_', group.list{group.idx}, '_', corr.list{corr.idx}, ...
        '_', p_value, '_', savetime);

    save(fullfile(path.outfolder, savename), 'results');
    saveas(gcf,fullfile(path.outfolder, ['sparse_disconns_' savename '.png']));

    %% 
    
    clearvars -except i_file infolder n_files path wkdir atlas group filename corr
    
end %i_file
