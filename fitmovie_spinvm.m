function fitmovie_spinvm(subject,arr,model,roiname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is for setting up parameters for model fitting with SPIN-VM
% arr = 1 for regularization parameter optimization
% arr = 2 for generating model weights and prediction scores
% model = 'sem' for category model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.lag = (-4:1:-2); % Hemodynamic delay
options.ratSamp = 1; % Ratio of training samples. e.g., 1/3 implies taking first 1/3*3600=1200 samples
options.tslots = 1;
options.regtype = 'L2'; % 'L2' for L2 regularization, 'L1' for L1 regularization
files.mask_size = 3;
files.cube_size = 3;
files.filter = ''; % Use '' for gaussian filter, 'avg' for average filter, 'log' for LoG filter
files.model = model;
fbase = './data/';
files.roiname = roiname;
files.roi = sprintf('%sROI_%s.mat',fbase,subject);
fout = sprintf('./results/spinvm_%s%s',lower(files.roiname),model); % Specify location for saving results

% Load stimulus matrix
files.stim = sprintf('%ssemgab2_movie_stim.mat',fbase);
lambda_first = 4; % This is related to the smallest regularization parameter (smallest exponent - 1)

if (arr == 1)
    options.cvfolds = 10; % Number of cross-validation folds
    options.vholdout = 10; % Number of chunks to hold out
    options.ctfolds = 10;
    options.tholdout = 10;
    % Set range of regularization parameters to use
    options.as = (2.^(lambda_first+1:lambda_first+10));    
    files.out  = sprintf('%s%s_Rv%s_%dc%d%s_mot.mat',fout,subject,files.roiname,files.cube_size,files.mask_size,files.filter);
    files.resp = sprintf('%ssem%s_presp_R.mat',fbase,subject); % BOLD response file
    disp(files.out); % Display filenames to make sure they are set correctly
    disp(files.stim);
    spinvm(files,1,options,subject);
    
elseif (arr > 1)
    options.cvfolds = 1000; % Number of cross-validation folds   
    % Below line would need to be changed for different set/range of regularization parameters
    lambda_constant = 1; % Regularization parameter multiplier
    files.reg = sprintf('./results/spinvm_%s%s%s_Rv%s_%dc%d%s_mot_cv_lambda.mat',lower(files.roiname),model,subject,files.roiname,files.cube_size,files.mask_size,files.filter);
    eval(['load ' files.reg]); % Load regularization parameters
    ts = squeeze(nanmean(nanmean(tccsall,2),1));
    numvx = size(ts,3);
    lambda_lap = zeros(1,numvx);
    lambda_reg = zeros(1,numvx);
    
    for ee=1:numvx
        tss = ts(:,:,ee);
        [~, indi] = max(tss(:));
        [a1, a3] = ind2sub(size(tss),indi);
        lambda_lap(ee) = a1; % Regularization parameters across neighborhoods (nei) optimizing prediction scores
        lambda_reg(ee) = a3; % Regularization parameters across features (feat) optimizing prediction scores
    end
    
    options.as_lap = lambda_constant*(2.^(lambda_lap+lambda_first));
    options.as_reg = lambda_constant*(2.^(lambda_reg+lambda_first));
    files.out  = sprintf('%s%s_Rv%s_%dc%d%s_mot.mat',fout,subject,files.roiname,files.cube_size,files.mask_size,files.filter);
    files.resp = sprintf('%ssem%s_presp_R.mat',fbase,subject);
    disp(files.out); % Display filenames to make sure they are set correctly
    disp(files.stim);
    spinvm(files,arr,options,subject);
end