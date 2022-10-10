%% ********************** POST-LQT SANITY CHECK ************************ %%
% 
% Written by Tamara Keßler, 07/2022
%
%%

% This script checks if the LQT ran correctly. This should best be checked
% for a subset of the data, as this script is not optimised and can take
% quite a while to run (sorry!)

%% Set up

clear; clc;

% Set up paths
% 1: Path to the older which contains all parcel_disconnection mat files 
% you want to check
path.infolder = 'D:\Tamara\DisconnectionMaps\Parcel_Disconnections\all\renamed';
folder = dir([path.infolder, '\*.mat']);
n_pat = size(folder,1);

% 2: Path to the atlas used for the parcellation in the LQT
path.atlas = 'C:\Users\Sperber\Desktop\BN\BNA_subregions_edited.xlsx';
atlas = readtable(path.atlas);
atlas_dim = size(atlas,1)*2;

path.data = 'D:\Tamara\Data\Patientlists_Groups';
RHLM_numbers = importdata(fullfile(path.data, 'all.txt'));

% Allocate some memory
percent_left_intra = zeros(n_pat,1);
is_sym = zeros(n_pat,1);
correct_dim = zeros(n_pat,1);

%%

% For every patient
for i_pat = 1:n_pat
    
    fprintf('Checking patient %d/%d\n',i_pat,n_pat);
    
    % Load parcel-wise disconnectivity matrix of patient
    filename = folder(i_pat).name;
    temp = load(fullfile(path.infolder, filename));
    matrix = temp.pct_sdc_matrix;
    
    % Check if matrix is symmetric
    is_sym(i_pat) = issymmetric(matrix);
    
    % Check if matrix has the same dimensions as the atlas
    dim_mat = size(matrix,1);
    correct_dim(i_pat) = isequal(atlas_dim,dim_mat); clear dim_mat
    
    % Find all non-zero elements in the matrix
    sparse_matrix = sparse(matrix);
    [row, col] = find(sparse_matrix);
    
    % Allocate some memory
    roi_disconns = cell(length(col),2);

    % For all disconnections in that patient
    for i_disconn = 1:length(col)

        % Even numbers are RH regions, odd are LH
        hem1 = mod(row(i_disconn),2);
        hem2 = mod(col(i_disconn),2);

        % If the numbers don't have the same parity, it's an
        % interhemispheric connection; else, it's intrahemispheric
        if hem1~=hem2
            roi_disconns{i_disconn,1} = 'inter';
            roi_disconns{i_disconn,2} = '---';
        else
            roi_disconns{i_disconn,1} = 'intra';
            
            % Check if even or odd number -> RH/LH (see above)
            if mod(hem1,2)==0
                roi_disconns{i_disconn,2} = 'right';
            else
                roi_disconns{i_disconn,2} = 'left';
                
            end %if-else (RH/LH)    
        end %if-else (inter- vs. intrahemispheric)
    end % i_disconn

    % Extract the hemispheres
    LR_intra = roi_disconns(:,2);
    
    counts_left = length(find(LR_intra=="left"));
    counts_right = length(find(LR_intra=="right"));
    counts_intra = length(find(roi_disconns(:,1)=="intra"));
    
    % Convert into categorical matrix for easier counting
%     cat_LR = categorical(LR_intra);
%     cat_LR_names = categories(cat_LR);
%     LR_counts = countcats(cat_LR);
%     
%     % Count how many left/right intrahemispheric connections there are
%     counts_left  = LR_counts(cat_LR_names=="left");
%     counts_right = LR_counts(cat_LR_names=="right");
%     counts_intra = counts_left + counts_right; clear counts_right
    
    % Calculate how many percent of the intrahemispheric connections are LH
    LR_percent = counts_left/counts_intra;
    
    left_intra_dcs = append(num2str(round(LR_percent*100)), '% left intra DCs (', ...
        num2str(counts_left), '/', num2str(counts_intra), ')');
    
    % to prevent errors
    if ~isempty(LR_percent)
        percent_left_intra(i_pat) = LR_percent;
    else
        percent_left_intra(i_pat) = 0;
    end
    
    clearvars -except i_pat folder n_pat path wkdir atlas ...
        percent_left_intra is_sym correct_dim atlas_dim RHLM_numbers
    
end %i_pat

%% Output

str1 = append('There are no left-intrahemispheric disconnections: ', ...
    num2str(~any(percent_left_intra)));
str2 = append('The matrices are symmetric: ', num2str(all(is_sym)));
str3 = append('The matrices have the same dimensions as the atlas: ', ...
    num2str(all(correct_dim)));

fprintf('\n');

fprintf('Results should all be "1"\n');
fprintf('%s\n',str1);
fprintf('%s\n',str2);
fprintf('%s\n',str3);

if any(percent_left_intra)
    idx = find(percent_left_intra);
    fprintf('There were %d left-intrahemispheric disconnections: \n',length(idx));

    
    for i_no = 1:length(idx)
        
        nr = idx(i_no);
        str4 = append('Patient ', RHLM_numbers{nr}, ' - ', ...
            num2str(round(percent_left_intra(nr)*100,2)), ' percent');
        
        fprintf('%s\n',str4);

    end % for
    
end % if