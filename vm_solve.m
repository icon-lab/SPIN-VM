function [ws, cmode] = vm_solve(X, Y, as, Xtest, testfeats)
% A simple multi-input multi-output ridge regressor
%    Input:
%          X: stimuli (SxN matrix, S is sample size, N is feature
%          size)
%          Y: responses (SxV matrix, V is model size)
%         as: an array of regularization parameters (Ax1 vector)
%      Xtest: testing data (), if model weights are not required
%  testfeats: the model features (X) that are used for obtaining predictions
%    Output:
%        ws: estimated weights (NxVxA matrix) if no testing data is provided
%            estimated response (TxVxA matrix) if testing data is given


if exist('Xtest', 'var')
    if ~exist('testfeats','var')
        testfeats = 1:size(X,2);
    end
end
if size(X,2)>size(X,1)
    cmode = 1;
elseif size(X,1) == size(X,2) && isequal(X,X')
    fprintf('X is symmetric. Assuming that this is a covariance matrix.\n');
    cmode = 2;
else
    cmode = 0;
end

switch cmode
    case 1
        [U, S] = eig(X*X');
    case 0
        [U, S] = eig(X'*X);
    case 2
        [U, S] = eig(X);
end
ds = diag(S);

if cmode
    U1 = X'*U;
    U2 = U'*Y;
else
    Uxy = U'*X'*Y;
end

if length(as) == size(Y,2)    
    % Do segmented ridge: different reg param for each voxel
    if ~exist('Xtest', 'var')
        ws = zeros(size(X,2),size(Y,2),1,'single');
    else
        ws = zeros(size(Xtest,1),size(Y,2),1,'single');
        if cmode
            U1 = (Xtest(:,testfeats)*U1(testfeats,:));
        else
            U =  (Xtest(:,testfeats)*U(testfeats,:));
        end
    end    
    uas = unique(as);
    for ii = 1:length(uas)
        % compute the models for a subset of voxels with the given optimal reg param.
        yind = find(as==uas(ii));        
        Sd = diag(1./(ds+uas(ii)));
        if cmode
            ws(:,yind,1) = U1*Sd*U2(:,yind);
        else
            ws(:,yind,1) = U*Sd*Uxy(:,yind);
        end
    end
else    
    % Do regular ridge: same param for all voxels
    if ~exist('Xtest', 'var')
        ws = zeros(size(X,2),size(Y,2),length(as),'single');
    else
        ws = zeros(size(Xtest,1),size(Y,2),length(as),'single');
        if cmode
            U1 = (Xtest(:,testfeats)*U1(testfeats,:));
        else
            U = (Xtest(:,testfeats)*U(testfeats,:));
        end
    end
    
    for ii = 1:length(as)        
        Sd = diag(1./(ds+as(ii)));        
        if cmode
            ws(:,:,ii) = U1*Sd*U2;
        else
            ws(:,:,ii) = U*Sd*Uxy;
        end
    end    
end
return