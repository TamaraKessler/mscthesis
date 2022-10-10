
%% Load binary sex
binsex = readmatrix('D:\Tamara\Data\sexbinary.csv');

%% Load disconnections

folder = 'D:\Tamara\Diskonnektionen\Parcel_Disconnections\all\renamed\parcels';

fprintf('Reading data...\n');
fileList = struct2cell(dir(fullfile(folder, '*.mat')));
temp=load(strcat(fileList{2,1}, '\', fileList{1,1})).pct_sdc_matrix;


images_2d = zeros(length(temp),length(temp),length(fileList)); %create blank for 2d images

for i=1:length(fileList)
    images_2d(:,:,i)=load(strcat(fileList{2,i}, '\', fileList{1,i})).pct_sdc_matrix;
    % IMPORTANT: same thing as in the previous comment!
end

%% create a 2d mask for the imaging data
mask=zeros(length(temp));
% the matrix is a symmetric matrix, so first remove all elements in the
% diagonal and below
for i=1:length(temp)
    for j=1:length(temp)
        if i<j
            mask(i,j)=1;
        end
    end
end

% then we need to remove rarely/never affected connections
images_2d_binary=double(images_2d>0);% binarise the image, every disconnection >0 is set to 1
image_2d_sum=sum(images_2d_binary,3);
image_2d_sum=image_2d_sum.*mask;

% get an overview on the number of patients with a disconnection for each
% ROI to ROI connection
% temp=sort(image_2d_sum(:),'descend');
% in 'temp', you can check how many connections are affected in at least x patients

% you can choose to only test connections that are affected in at least threshold = x
% patients (I set it to 5, but feel free to change it)
threshold = 5;
mask(image_2d_sum<threshold)=0;
mask_vect=mask(:);

% vectorise the 2d images
images_vect=zeros(length(fileList),sum(mask(:))); %create blank
for i=1:length(fileList)
    temp=images_2d(:,:,i);
    temp=temp(:);
    images_vect(i,:)=temp(mask_vect==1);
end

%% split into F and M

images_vect_F = zeros(length(fileList),sum(mask(:))); 
images_vect_M = zeros(length(fileList),sum(mask(:)));

for i_pnum = 1:length(fileList)
   
    if binsex(i_pnum)==1
        images_vect_F(i_pnum,:) = images_vect(i_pnum,:);
    elseif binsex(i_pnum)==0
        images_vect_M(i_pnum,:) = images_vect(i_pnum,:);
    else
        continue
    end
    
end

% delete rows with all zeros
images_vect_F(~any(images_vect_F,2),:) = [];
images_vect_M(~any(images_vect_M,2),:) = [];

%%

% allocate some memory
p_orig = zeros(sum(mask(:)),1);
h_orig = zeros(sum(mask(:)),1);
%stats_orig=zeros(sum(mask(:)),1);

% for every disconnection
for i_discon = 1:size(p_orig,1)

    % use wilcoxon test (F vs M)
    [p,h] = ranksum(images_vect_F(:,i_discon),images_vect_M(:,i_discon));
    p_orig(i_discon) = p;
    h_orig(i_discon) = h;
    %stats_orig(i_discon) = stats;

end

%%

num_sig = length(find(h_orig==1));