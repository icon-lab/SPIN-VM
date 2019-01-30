%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SPIN-VM DEMO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script demonstrates the fMRI modeling technique developed in the
% following paper called spatially informed voxelwise modeling (SPIN-VM)
% Çelik, E., Dar, S. U. H., Yılmaz, Ö., Keleş, Ü., & Çukur, T. (2019). 
% Spatially informed voxelwise modeling for naturalistic fMRI experiments. NeuroImage, 186, 741-757.
% doi: 10.1016/j.neuroimage.2018.11.044
% SPIN-VM leverages correlations across neighboring voxels to better
% predict fMRI responses obtained on independent test data and to better
% capture locally congruent information representations across cortex compared to VM.

% Initialize variables
sub = 'S1'; % subject name
roiname = 'PPA'; 
cube_size = 3;
mask_size = 3;
filter_type = 'gaussian';
model = 'sem'; % 'sem' for category model, 'gab314' for motion-energy model
fdata = './data/';
mkdir results % Create a folder for storing results
fresults = './results/';

% Load ROI
froi = sprintf('%sROI_%s.mat',fdata,sub); 
eval(['load ' froi]);
roiloc = strfind(roilist,roiname);
roiid = find(not(cellfun('isempty',roiloc)));
targetvoxels = roivox{roiid}; % Indices of the voxels in the specified ROI

% Generate a neighborhood matrix and a mask matrix that stores weights to be assigned to the neighbors
neighbor_searchlight(sub,cube_size,mask_size,filter_type);

% Generate a Laplacian matrix using the neighborhood and mask matrices obtained in the previous step
roi_neighbors_wb(sub,roiname,cube_size,mask_size);

% Compute Schur decomposition of the Laplacian matrix and save it for further use
schur_precompute(sub,roiname,cube_size,mask_size);

% Optimize regularization parameters for VM for the subject/model/ROI specified
fitmovie_vm(sub,1,model,roiname);
% Generate model weights and prediction scores using the optimal regularization parameters for VM
fitmovie_vm(sub,2,model,roiname);

% Optimize regularization parameters for SPIN-VM for the subject/model/ROI specified
fitmovie_spinvm(sub,1,model,roiname);
% Generate model weights and prediction scores using the optimal regularization parameters for SPIN-VM
fitmovie_spinvm(sub,2,model,roiname);

% Load prediction scores for both VM and SPIN-VM
tccs_ind_vm = h5read(sprintf('%s%s%s%s_Rv%s_matched_motcvall_corrs.hf5',fresults,lower(roiname),model,sub,roiname),'/tccs_ind');
tccs_ind_spinvm = h5read(sprintf('%sspinvm_%s%s%s_Rv%s_%dc%d_motcvall_corrs.hf5',fresults,lower(roiname),model,sub,roiname,cube_size,mask_size),'/tccs_ind');

% Load correction factors
% Raw correlation coefficients are biased downward by noise in the measured
% BOLD responses. Therefore, they were corrected for noise using these
% correction factors. See below paper for details.
% David, S. V., & Gallant, J. L. (2005). 
% Predicting neuronal responses during natural vision. Network: Computation in Neural Systems, 16(2-3), 239-260.
% doi: 10.1080/09548980500464030
cfact = h5read(sprintf('%ssem%s_Rvallcvall_noise_cf.hf5',fdata,sub),'/cfact');
cfact = cfact';
cfact = cfact(targetvoxels);
cfact(isnan(cfact))=1; % Make the NaN values in cfact 1

scores_spinvm = transpose(mean(tccs_ind_spinvm,2)); % Take the mean to obtain a single prediction score value for each voxel
scores_vm = transpose(mean(tccs_ind_vm,2));
scores_spinvm = scores_spinvm.*cfact; % Applying correction factors to prediction scores
scores_vm = scores_vm.*cfact;
scores_spinvm(scores_spinvm>1) = 1; % Set prediction scores to 1 if they are higher than 1 due to cfacts
scores_vm(scores_vm>1) = 1; 
scores_spinvm(scores_spinvm<0) = 0; % Set prediction scores to 0 if they are lower than 0
scores_vm(scores_vm<0) = 0;

% Plot prediction scores for VM and SPIN-VM for comparison
xx = 0:0.001:1;
figure
plot(scores_vm,scores_spinvm,'o')
hold on
plot(xx,xx,'r')
title('Prediction Scores')
xlabel('VM')
ylabel('SPIN-VM')
