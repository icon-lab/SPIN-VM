function roi_neighbors_wb(sub,roiname,cube_size,mask_size)
% This script generates and saves a graph Laplacian matrix
% Example => roi_neighbors_wb('S1','PPA',3,3)
fbase = './data/';

% load ROI mask

rois = sprintf('%sROI_%s.mat',fbase,sub);
eval(['load ' rois]); % load ROIs
roiloc = strfind(roilist,roiname);
roiid = find(not(cellfun('isempty',roiloc)));
targetvoxels = roivox{roiid}; % targetvoxels stores the indices of the voxels in the specified ROI

fbase = './data/'; % path where neighborhood file is located
load(sprintf('%sneighbor_%s_%dc%d.mat',fbase,sub,cube_size,mask_size),'position_matrix','mask_matrix')
nei = position_matrix(targetvoxels,:);
mask = mask_matrix(targetvoxels,:);
[~,b] = ismember(nei,targetvoxels);
sz = size(nei);
c_nei = single(zeros(sz(1))); % c_nei is denoted as matrix 'C' in the paper
T = single(zeros(sz(1)));

for vox = 1:sz(1)
    [~,d] = find(b(vox,:)~=0);
    c_nei(vox,b(vox,d)) = mask(vox,d);
end
c_nei(find(eye(length(T))==1)) = 0;
T(find(eye(length(T))==1)) = sum(c_nei);
L = T - c_nei; % graph Laplacian

% Only need to save L
eval(['save -v7.3 ' sprintf('%s_roi_neighbor_%s_%s_%dc%d',fbase,sub,roiname,cube_size,mask_size)  '.mat L']);
