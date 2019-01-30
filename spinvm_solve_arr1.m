function [Y_pred] = spinvm_solve_arr1 (X,Y,U,S,reg,X_test)
% This script computes and outputs predicted responses (Y_pred) for each pair of reg. param.
% It first finds associated model weights (w) by solving AW+WB = M, where
% A = auto-covariance matrix of size (3 x Nfeat)x(3 x Nfeat)
% B = Laplacian matrix of size (Nvox)x(Nvox), where B = lambda_nei*L
% M = cross-covariance matrix of size (3 x Nfeat)x(Nvox)
% See 'Pseudocode for SPIN-VM' from the paper for more details
X_test = single(X_test);
M = -X'*Y;
Mtilde = M*U;
Art = X'*X;
IA = speye(length(Art));
Y_pred = single(zeros(size(X_test,1),size(U,1),length(reg),length(reg)));
for i = 1:length(reg) % lambda_feat
    A = (Art+reg(i)*IA); % X'*X+lambda*I
    [Q, D] = eig(A); % Find eigenvalues and eigenvectors of X'*X+lambda*I
    D = diag(D);
    Qt = Q';
    sizA = size(A,1);
    Drep = repmat(D,1,size(U,1));
    for j = 1:length(reg) % lambda_nei
        P = spinvm_rep(reg(j)*S,Drep,sizA);
        Xtilde = -Q*(P.*(Qt*Mtilde));
        w = Xtilde*U';
        w = real(w);
        w = single(w); % Save as single to reduce size
        Y_pred(:,:,j,i) = X_test*w; % Calculate predicted response
    end
end
