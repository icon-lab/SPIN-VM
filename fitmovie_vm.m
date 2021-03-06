function fitmovie_vm(subject,arr,model,roiname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is for setting up parameters for model fitting with VM
% arr = 1 for regularization parameter optimization
% arr = 2 for generating model weights and prediction scores
% model = 'sem' for category model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.lag = (-4:1:-2); % Hemodynamic delay
options.ratSamp = 1; % Ratio of training samples. 1/3 implies taking first 1/3*3600=1200 samples
options.tslots = 1;
options.regtype = 'L2'; % 'L2' for L2 regularization, 'L1' for L1 regularization
files.model = model;
fbase = './data/';
files.roiname = roiname;
files.roi = sprintf('%sROI_%s.mat',fbase,subject);
fout = sprintf('./results/%s%s',lower(files.roiname),model); % Specify location for saving results

% Load stimulus matrix
files.stim = sprintf('%ssemgab2_movie_stim.mat',fbase);
lambda_first = 4; % This is related to the smallest regularization parameter (smallest exponent - 1)

if (arr == 1)
    options.cvfolds = 10; % Number of cross-validation folds
    options.vholdout = 10; % Number of chunks to hold out
    options.ctfolds = 10;
    options.tholdout = 10;
    % Set range of regularization parameters to use
    options.as = 2.^(lambda_first+1:lambda_first+10);      
    files.out  = sprintf('%s%s_Rv%s_matched_mot.mat',fout,subject,files.roiname);      
    files.resp = sprintf('%ssem%s_presp_R.mat',fbase,subject); % BOLD response file
    disp(files.out);
    disp(files.stim);
    vm(files,1,options);
    
elseif (arr > 1)    
    options.cvfolds = 1000; % Number of cross-validation folds
    files.reg = sprintf('./results/%s%s%s_Rv%s_matched_mot_cv_lambda.mat',lower(files.roiname),model,subject,files.roiname);
    eval(['load ' files.reg]);
    t1 = squeeze(nanmean(nanmean(tccsall,2),1));
    [~, mi1] = max(t1,[],1);
    options.as = as(mi1);
    files.out  = sprintf('%s%s_Rv%s_matched_mot.mat',fout,subject,files.roiname);
    files.resp = sprintf('%ssem%s_presp_R.mat',fbase,subject);
    disp(files.out);
    disp(files.stim);
    vm(files,arr,options);
end