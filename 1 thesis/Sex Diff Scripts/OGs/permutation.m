function [data_perm, idx_perm] = permutation(data, stream)
% ***************************** Permutation *******************************
%
% Lisa Röhrig, 18-12-2020, last update on 09-07-2021
%
% Permutation process of input data
% Shuffles randomly data

idx_data = (1:size(data, 1))';
idx_perm = randperm(stream, length(idx_data))';
data_perm = data(idx_perm, :);

end