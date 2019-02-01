function [position_matrix,mask_matrix]=neighbor_searchlight(sub,cube_length,mask_length,filter_type)
% Input: sub = name of subject, cube_length and mask_length = window size (optimally both the same and equal to 3) 
% filter_type = 'gaussian', 'average', or 'log'
% Output: a neighborhood (position) matrix and a mask matrix that stores weights to be assigned to the neighbors
% Example => neighbor_searchlight('S1',3,3,'gaussian')

fbase = './data/';
load(sprintf('%ssem%s_presp_R.mat',fbase,sub),'tvoxels'); % tvoxels = array storing the indices of brain voxels
slices = 32; % number of transverse slices

offset = 2*cube_length;
temp_padded = zeros(100+offset,100+offset,slices+offset);
temp_not_padded = zeros(100,100,slices);
temp_not_padded(tvoxels) = tvoxels;
temp_padded((1:100)+cube_length,(1:100)+cube_length,(1:slices)+cube_length) = temp_not_padded;
window = (cube_length-1)/2;
neighborhood = ((cube_length^3))+1;
position_matrix = zeros(length(tvoxels),neighborhood);
mask_matrix = zeros(length(tvoxels),neighborhood);
mask = fspecial3(filter_type, [mask_length mask_length mask_length]);
mask_center = (mask_length+1)/2;
mask = mask(mask_center-window:mask_center+window,mask_center-window:mask_center+window,mask_center-window:mask_center+window);
temp_matrix = zeros(size(mask));

for voxel_index=1:length(tvoxels)    
    
    [row,column,slice] = ind2sub([100,100,slices],tvoxels(voxel_index)); % each slice is of size 100x100       
    % check the range of neighbors    
    s_l = slice-window+cube_length;
    s_r = slice+window+cube_length;   
    x_l = row-window+cube_length;
    x_r = row+window+cube_length;    
    y_l = column-window+cube_length;
    y_r = column+window+cube_length;    
    % generate a temporary neighborhood matrix
    temp_matrix = temp_padded(x_l:x_r,y_l:y_r,s_l:s_r);
    [~,~,non_zero] = find(temp_matrix);
    % apply mask
    non_zero_mask = find(temp_matrix~=0);
    p = zeros(cube_length,cube_length,cube_length);
    p(non_zero_mask) = 1;
    non_zero_mask = p.*mask;
    [~,~,non_zero_mask] = find(non_zero_mask);   
    position_matrix(voxel_index,1:length(non_zero)) = non_zero';
    mask_matrix(voxel_index,1:length(non_zero)) = non_zero_mask';
end
[~,position_matrix] = ismember(position_matrix,tvoxels);

if strcmp(filter_type,'average')
    save(sprintf('%sneighbor_%s_avg%dc%d.mat',fbase,sub,cube_length,mask_length),'position_matrix','mask_matrix','-v7.3')
elseif strcmp(filter_type,'gaussian')
    save(sprintf('%sneighbor_%s_%dc%d.mat',fbase,sub,cube_length,mask_length),'position_matrix','mask_matrix','-v7.3')
elseif strcmp(filter_type,'log')
    save(sprintf('%sneighbor_%s_log%dc%d.mat',fbase,sub,cube_length,mask_length),'position_matrix','mask_matrix','-v7.3')
end