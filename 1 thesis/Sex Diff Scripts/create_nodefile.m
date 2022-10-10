temp = load('D:\Lesion_Quantification_Toolkit\Support_Tools\Parcellations\Schaefer_Yeo\Plus_Subcort\Schaefer2018_100Parcels_7Networks_order_plus_subcort.mat');
yeo_100_7 = temp.t;

n_nodes = size(yeo_100_7,1);

nodes = zeros(n_nodes,5);
nodes(:,1) = yeo_100_7.X;
nodes(:,2) = yeo_100_7.Y;
nodes(:,3) = yeo_100_7.Z;
nodes(:,4) = repelem(1,n_nodes);
nodes(:,5) = repelem(1,n_nodes);

save newnodes.txt nodes -ascii -double