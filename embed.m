function [eX, timeidx, tauidx] = embed(X,tau,fskip)
% This script is used for stimulus preprocessing
% e.g., incorporating hemodynamic delays and skipping time frames
if (nargin<3)
    fskip.flag = 0;
    fskip.mode = 2;
end

if isfield(fskip,'mode')
    mode = fskip.mode;
else
    mode = 2;
end

% embed the first signal in its temporal context
[D, T] = size(X);

if mode == 2
    startInd = abs(tau(1)) + 1;
    stopInd	= T - abs(tau(end));
    len	= stopInd - startInd + 1;
elseif mode == 1
    startInd = 1;
    stopInd	= T;
    len	= stopInd - startInd + 1;
end

if (mode == 1)
    nlag = length(tau);
    eX = zeros(nlag,D,T);
    for shind = 1:nlag
        eX(shind,:,:) = circshift(X,[0 -tau(shind)]);
    end
    eX = reshape(eX,nlag*D,T);
    eX = eX(:,startInd:stopInd);
end

if (mode == 2)
    % create a column vector that contains the indices of the first segment
    idx = repmat((startInd:stopInd)', 1, length(tau)) + repmat(tau, len, 1);
    % create (linear) indices for the different dimensions
    dim_offset = repmat( (0:D-1)*T, length(tau)*len, 1);
    idx = repmat(idx(:), 1, D) + dim_offset;
    % for the linear indices we need column-signals
    X = X';
    % get the data (D channels, segments are concatenated) and reshape it
    eX = reshape(X(idx), len, length(tau)*D)';
end
tauidx = repmat(tau',D,1);

if (fskip.flag == 1)
    %skip nskip time frames at the beginning of each run
    TR = fskip.TR;
    nskip = 6/TR;
    if (fskip.val == 1)
        frtime = 60/TR;
    else
        frtime = 600/TR;
    end
    baseslot = (nskip+1):1:frtime;
    idSamples = [];
    for nslot = 1:ceil(T/frtime)
        idSamples = [idSamples (baseslot + frtime*(nslot-1))];
    end
    subInd = find((idSamples>=startInd) & (idSamples<=stopInd));
    skipInd = idSamples(subInd);
    eX = eX(:,skipInd-startInd+1);
    timeidx = skipInd;
else
    timeidx = startInd:stopInd;
end

