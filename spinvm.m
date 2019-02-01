function spinvm(files, stage, options, sub)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spatially informed voxelwise modeling for predefined ROIs or whole brain
% --
% files         -  names of output, stim, resp and testcc files
% stage         -  stage of processing
%                  (1 = Reg. Param Selection)
%                  (2 = Model Estimation)
% options       -  analysis parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------
% Parse inputs
%-------------------------------------------------------
if (nargin < 1)
    error('Need to input a filename@ Quitting...');
elseif (nargin < 2)
    error('Need to indicate processing stage! Quitting...');
elseif (nargin < 3)
    error('Options missing from input! Quitting...');
end

%-------------------------------------------------------
% Ready to load data
%-------------------------------------------------------
% Load preprocessed stimuli
eval(['load ' files.stim]);

% Apply PCA on stimulus matrix
[~, spca, ~, ~, exx] = pca(stim(:,2:end));
npca = 300; % number of principal components used
spca = spca(1:3870,1:npca);
sum(exx(1:npca)) % explained by first npca components
stim = [stim(:,1),spca]; % Append motion energy column

% Load voxel responses
eval(['load ' files.resp]);

%-------------------------------------------------------
% Set default options
%-------------------------------------------------------
if (isfield(options,'regtype'))
    regtype = options.regtype;
else
    regtype = 'L2';
end
if (isfield(options,'cvfolds'))
    cvfolds = options.cvfolds;
else
    cvfolds = 1;
end
if (isfield(options,'vholdout'))
    vholdout = options.vholdout;
else
    vholdout = 5; 
end
if (isfield(options,'ctfolds'))
    ctfolds = options.ctfolds;
else
    ctfolds = 1;
end
if (isfield(options,'tholdout'))
    tholdout = options.tholdout;
else
    tholdout = 3;
end
if (isfield(options,'ratSamp'))
    ratSamp  = options.ratSamp;
else
    ratSamp  = 1;
end
if (isfield(options,'trnSamp'))
    trnSamp  = options.trnSamp;
else    
    trnSamp = 1:length(ranges.crange);    
end
if (isfield(options,'lag'))
    lag  = options.lag;
else
    lag  = 0;
end
if (isfield(options,'as'))
    as      = options.as;
else
    as      = 100*(2.^(0:15)); % default reg. parameters
end

%------------------------
% Voxel selection (ROI)
%------------------------

% load mask-voxel information
disp('Running on ROI...');
eval(['load ' files.roi]);
roiloc = strfind(roilist,files.roiname);
roiid = find(not(cellfun('isempty',roiloc)));
targetvoxels = roivox{roiid}; % Indices of the voxels in the specified ROI
numVox = length(targetvoxels);

fprintf('Running on %d voxels...\n',numVox);
%-------------------------------------------------------
% prepare the training stimulus and response with lag
%-------------------------------------------------------
if strcmp(files.model,'sem') && stage==2
    stim = stim(:,2:end); % Remove motion-energy feature for model estimation
end
stim = single(stim);
numFeat  = size(stim,2); % Number of features
% Clean NaN entries in the response data
resp(find(isnan(resp))) = 0;
resp = single(resp);
% Restrict time frames
numSamples = floor(ratSamp*length(trnSamp));
numTrainingSamples = (1:numSamples);
% Frame skipping
fskip.TR   = 2; % TR for acquisition (2 s)
fskip.flag = 0; 
fskip.val  = 0; 
% Stimulus prep and demeaning
fskip.mode = 1; 
[stimsh, timeidx, ~] = embed(stim(trnSamp(numTrainingSamples),:).',lag,fskip);
stimsh = stimsh.';
stimsh = stimsh - repmat(mean(stimsh,1),size(stimsh,1),1);
% Response prep and demeaning
respsh = resp;
respsh = respsh(timeidx,:);
respsh = respsh - repmat(mean(respsh,1),size(respsh,1),1);

% For cross-validation, seperate the data into cvfolds chunks
numTrainingSamples = 1:length(timeidx);
if stage == 1
    % find the start-end points of the chunks
    chunklen = 25; % 1 chunk = 2 s, 25 chunks = 50 s
    numchunk = round(length(numTrainingSamples)/chunklen);
    mid_frames = [round((0:(numchunk-1))*length(numTrainingSamples)/numchunk)+1 ...
        length(numTrainingSamples)+1];
    % divide the data into test and training sets given cvfolds
    sblock = 1:numchunk;
    numvblock = vholdout;
    for cvidx = 1:cvfolds
        if (length(sblock)<numvblock), sblock = 1:numchunk; end % restart if sblock is done
        testblks = randsample(sblock,numvblock);
        trnblks{cvidx} = setdiff(1:numchunk,testblks);
        sblock = setdiff(sblock,testblks);
        tsample = [];
        for blkidx = 1:numvblock
            tsample = cat(2,tsample,mid_frames(testblks(blkidx)):(mid_frames(testblks(blkidx)+1)-1));
        end
        testingSamples{cvidx}  = tsample;
        trainingSamples{cvidx} = setdiff(numTrainingSamples, testingSamples{cvidx});
    end
    
elseif stage == 2 
    %-------------------------------------------------------
    % prepare the validation stimulus and response with lag
    %-------------------------------------------------------
    fskip.flag = 0; % Turn on frame skipping
    fskip.val  = 1; % Set to 1 for validation data
    [stimshv, timeidx, ~] = embed(stim(ranges.vrange,:).',lag,fskip);
    stimshv = stimshv.';
    stimshv = stimshv - repmat(mean(stimshv,1),size(stimshv,1),1);
    % Prep response matrix     
    respshv = resp(ranges.vrange,:);    
    respshv = respshv(timeidx,:);
    respshv = respshv - repmat(mean(respshv,1),size(respshv,1),1);
    % Extract validation samples
    valSamples = 1:length(timeidx);
end

% Load Schur decomposition variables
if strcmp(files.filter,'avg')
    load(sprintf('./data/%sschurvarsavg%s%dc%d.mat',sub, files.roiname, files.cube_size, files.mask_size),'U','S')
elseif strcmp(files.filter,'log')
    load(sprintf('./data/%sschurvarslog%s%dc%d.mat',sub, files.roiname, files.cube_size, files.mask_size),'U','S')
else
    load(sprintf('./data/%sschurvars%s%dc%d.mat',sub, files.roiname, files.cube_size, files.mask_size),'U','S')
end

switch stage
    
    %-------------------------------------------------------
    % Run cross-validation for determining the reg. param
    %-------------------------------------------------------
    case {1}        
        % initialize parameters
        tccsall = zeros(cvfolds,ctfolds,length(as),length(as),length(targetvoxels));
        mccsall = zeros(cvfolds,ctfolds,length(as),length(as));
        
        % Generate prediction scores across cross-validation folds
        for cvidx = 1:cvfolds
            fprintf('iter no %d of %d(cv)\n',cvidx,cvfolds);
            %Get the training blocks for that cvfold and resample from them
            sblock = trnblks{cvidx};
            numtblock = tholdout;
            for tidx = 1:ctfolds
                fprintf('%d of %d(ct)\n',tidx,ctfolds);
                if (length(sblock)<numtblock), sblock = trnblks{cvidx}; end % restart if sblock is done
                testblks = randsample(sblock,numtblock);
                sblock = setdiff(sblock,testblks);
                tsample = [];
                for blkidx = 1:numtblock
                    tsample = cat(2,tsample,mid_frames(testblks(blkidx)):(mid_frames(testblks(blkidx)+1)-1));
                end
                subtstSamples{tidx} = tsample;
                subtrnSamples{tidx} = setdiff(trainingSamples{cvidx}, subtstSamples{tidx});
                
                % Extract testing and training data
                Xtest = stimsh(subtstSamples{tidx},:);                
                Ytest = respsh(subtstSamples{tidx},:);                
                Xtrain = stimsh(subtrnSamples{tidx},:);
                Ytrain = respsh(subtrnSamples{tidx},:);
                
                % Obtain predicted responses
                Ytest_pred = spinvm_solve_arr1(Xtrain, Ytrain, U, S, as, Xtest);
                Ytest_pred = Ytest_pred(:,:,:,:);                
                
                % Calculate prediction scores
                for asidx1 = 1:length(as)
                    for asidx2 = 1:length(as)
                        % Calculate corr for all voxels
                        tccs = mvn_corr_fast(Ytest, Ytest_pred(:,:,asidx1,asidx2));
                        tccsall(cvidx,tidx,asidx1,asidx2,:) = tccs;
                        mccs = nanmean(tccs);
                        % Store all correlation values
                        mccsall(cvidx,tidx,asidx1,asidx2) = mccs;
                    end
                end
            end
        end
        
        % Save the output as single to reduce footprint
        tccsall = single(tccsall);
        mccsall = single(mccsall);
        eval(['save -v7.3 ' files.out(1:end-4) '_cv_lambda.mat as tccsall  mccsall trainingSamples']);
        %-----------------------------------------------------------------
        % Plot prediction scores across reg. params to see optimal values
        %-----------------------------------------------------------------
        mccsall = reshape(mccsall, [cvfolds*ctfolds length(as) length(as)]);
        
        [cm, ci] = max(mean(mccsall,1));
        fprintf('Best reg. param based on corrs: %d\n', as(ci));
        
        mcc = mean(mccsall,1);
        mcc = reshape(mcc,[length(as) length(as)]);
        
        figure
        semilogx(as, mean(mcc,1));
        xlabel('regularization parameter');
        ylabel('average correlation coefficient');
        title(sprintf('Max mean test cc: %.3f', cm));

    case {2}
        % Number of time lags
        nTau = length(lag);
        tccs_ind = zeros(numVox,cvfolds);
        % select the indices of the test features used in prediction        
        sfeats = 1:numFeat;
        testfeats = (nTau*(sfeats(1)-1)+1):(nTau*sfeats(end));
                
        % first fit the model to all training data with optimal reg. param.
        
        % training data
        Xtrain = stimsh;
        Ytrain = respsh;
        % validation data
        Xvalidation = stimshv(valSamples,:);
        Yvalidation = respshv(valSamples,:);
        clear resp respsh stim stimsh
        % Obtain model weights using optimal reg. parameters
        ws = spinvm_solve_arr2(Xtrain, Ytrain, U, S, options.as_lap, options.as_reg);        
        
        % Obtain predicted responses for independent validation data
        Yval_pred = Xvalidation(:,testfeats)*ws(testfeats,:);
        
        % lag-averaged voxel strfs
        vsrf = reshape(ws, [nTau, numFeat, numVox]);
        vsrf = single(squeeze(mean(vsrf,1)));
        vsrf = vsrf(:,1:numVox);
        % save the model
        hdf5write(sprintf('%s_ridge_srfs.hf5',files.out(1:end-4)),'vsrf',vsrf,'ws',ws);
        clear vsrf;
                
        % Obtain prediction scores across cross-validation folds
        totaldur = length(valSamples);
        for cvidx = 1:cvfolds            
            timeind = randsample(totaldur,round(totaldur*0.8),false); % with replacement
            tccs_ind(:,cvidx) = mvn_corr_fast(Yvalidation(timeind,:), Yval_pred(timeind,:));
            mccs_ind(cvidx) = nanmean(tccs_ind(:,cvidx));            
            fprintf('avg. cc for ind: %.5f\n',mccs_ind(cvidx));            
        end        
        % Save all correlation values
        tccs_ind = single(tccs_ind);
        mccs_ind = single(mccs_ind);
        hdf5write(sprintf('%scvall_corrs.hf5',files.out(1:end-4)),'tccs_ind',tccs_ind,'mccs_ind',mccs_ind);
end