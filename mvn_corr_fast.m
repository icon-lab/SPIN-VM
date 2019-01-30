function predcc = mvn_corr_fast(X, Y)
% Compute pairwise correlations between identical columns of X and Y

if ~~sum(isnan(Y(:))|isnan(X(:)))
    error('NaN entries!')
end

predcc = sum((zscore(X)).*(zscore(Y)))/(size(X,1)-1);