
% Set working directory
wkdir = "D:\Tamara";

%% Get user input: Which task should be done - renaming or moving files?

prompt.task = "Select: 1 = Rename Files or 2 = Move Files";
dlg.task = "What do you want to do?";
dims = [1 60];
default = "1";

answer = inputdlg(prompt.task, dlg.task, dims, default);
task = str2num(answer{1});


%% Get user input: Select input folder & input file type

prompt.infolder = "Select input folder";
path.infolder = uigetdir(wkdir, prompt.infolder);

ft.list = {'.nii', '.mat'};
ft.idx = listdlg('PromptString','Select the input file type', ...
    'SelectionMode','single','ListString',ft.list);

infolder = dir([path.infolder, '\*', ft.list{ft.idx}]);
n_pat = size(infolder,1);

%% Get user input: Get patient list

prompt.groupfile = "In the next step, please select the .txt file containing relevant patient numbers";
f = msgbox(prompt.groupfile); pause;
[filename.patnum, path.patnum] = uigetfile('*.txt');

patnums = importdata(fullfile(path.patnum, filename.patnum));

if n_pat ~= size(patnums,1)
    error("Number of files doesn't match provided list of patients");
end

%% Perform selected task

% If task = 1 / rename files
if task == 1
    
    suffix.list = {'(None)', '_disconmap', '(Other)'};
    suffix.idx  = listdlg('PromptString','Select a suffix for the files', ...
        'SelectionMode','single','ListString',suffix.list);
    
    % For all patients
    for i_pat = 1:n_pat
        
        % Create new filename
        if suffix.idx == 1
            new_name = fullfile(path.infolder, [RHLM_numbers{i_pat}, ft.list{ft.idx}]);
        elseif suffix.idx == 2
            new_name = fullfile(path.infolder, [RHLM_numbers{i_pat}, ...
                suffix.list{suffix.idx}, ft.list{ft.idx}]);
        end %if-else
        
        % Retrieve old filename
        old_name = fullfile(path.infolder, infolder(i_pat).name);
        
        % Rename the file (by "moving" it within the same folder)
        movefile(old_name,new_name);
        
        clear old_name new_name
        
    end %i_pat
    
    % If task = 2 / move files
elseif task == 2
    
    % Get user input: Select output folder
    prompt.outfolder = "Select output folder";
    path.outfolder = uigetdir(wkdir, prompt.outfolder);
    
    % For all patients
    for i_pat = 1:n_pat
        
        % Retrieve filename & path
        pat_name = infolder(i_pat).name;
        old_path = fullfile(path.infolder, pat_name);
        
        % Move file to selected output folder
        movefile(old_path,path.outfolder);
        
        clear pat_name old_path
        
    end %i_pat
    
else
    
    error("Selected task type invalid");
    
end %if-else