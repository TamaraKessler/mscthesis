%% ****************  Check Bounding Box and Max Value ********************
%-----------------------------------------------------------------------
% Lisa Röhrig, 27-07-2020; adjusted by Tamara Keßler (06/2022)
%-----------------------------------------------------------------------
%
% Controls for same bounding boxes of image sample and correct data type
% (maximum value of image matrix has to be 1 instead of 255...).

clear; clc;

%filedir = 'D:\Lisa\Lesions - Lisa\';
filedir = 'D:\Tamara'
cd(filedir);

% ........................... Select Files ................................

% User selects one or multiple files
[files, path] = uigetfile('*.nii', 'Pick one or multiple lesion maps.', 'MultiSelect', 'on');

% If only one file is selected, transfer it into a cellarray
if ~iscell(files)
    files = {files};
end % Now files is a cell array regardless of the number of selected files.

% ............................ CONTROL CHECK ..............................

fprintf('Control Check (Bounding Box, Maximum Value, Minimum Value) .....\n');

% Create empty list
List_BB = [];
x = 1;
List_BB{x, 1} = 'Image_file';
List_BB{x, 2} = 'Bounding_Box';
List_BB{x, 3} = 'Max_Value';
List_BB{x, 4} = 'Min_Value';
x = x+1;

MaxVal_all = 0; MinVal_all = 0; BB_cmp = 1;

% For all selected files
for i = 1:length(files)
        
    % Image file
    image_header = spm_vol([path files{i}]);
    image_file = spm_read_vols(image_header);
    List_BB{x, 1} = files{i};

    % Bounding Box
    BB = size(image_file);
    List_BB{x, 2} = [num2str(BB(1)) 'x' num2str(BB(2)) 'x' num2str(BB(3))];
    
    % Maximum Value of Image Matrix
%    MVal_pre = max(image_file, [], 'all');
    
    % Change Maximum Value 255 to 1
%     if MVal_pre > 1
%         image_file(image_file > 1) = double(1);
%         fprintf('\nMaximum value of image matrix was greater than 1, now it is corrected to 1.');
%         fprintf('\nImage file:\t%s\n', files{i});
%         spm_write_vol(image_header, image_file);
%     end
    MaxVal = max(image_file, [], 'all');
    List_BB{x, 3} = MaxVal;
    MinVal = min(image_file, [], 'all');
    List_BB{x, 4} = MinVal;
    
    x = x+1;
end

% Comparison
for i = 2:(length(List_BB)-1)
    
    % Maximum value
    if List_BB{i, 3} > MaxVal_all
        MaxVal_all = List_BB{i, 3};
    end
    
    % Minimum value
    if List_BB{i, 4} < MinVal_all
        MinVal_all = List_BB{i, 4};
    end
    
    % Same Bounding Boxes?
    if BB_cmp ~= 0
        if strcmp(List_BB{i, 2}, List_BB{i+1, 2}) == 1
            BB_cmp = 1; % True
        else
            BB_cmp = 0; % False
            break
        end
    end
        
end

% Output
if BB_cmp == 1
    fprintf('\nAll selected image files have the same bounding boxes: \tTRUE');
    fprintf('\nBounding Box: \t%d x %d x %d\n', BB(1), BB(2), BB(3));
else
    fprintf('\nAll selected image files have the same bounding boxes: \tFALSE\n');
end
fprintf('\nThe maximum value of all image matrices is: \t%d\n', MaxVal_all);
fprintf('\nThe minimum value of all image matrices is: \t%d\n', MinVal_all);

% Print list?
fprintf('\nWant to print the whole list?\n');
a = [];
while 1
    a = input('Yes (1) or No (0): ');
    if a == 1
        display(string(List_BB));
        break
    elseif a == 0
        break
    end
end

