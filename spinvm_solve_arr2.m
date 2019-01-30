function [ws] = spinvm_solve_arr2 (X,Y,U,S,reg_lap,reg_reg)
% This script computes and outputs model weights (ws) for all voxels
% It first finds associated model weights (w) by solving AW+WB = M, where
% A = auto-covariance matrix of size (3 x Nfeat)x(3 x Nfeat)
% B = Laplacian matrix of size (Nvox)x(Nvox), where B = lambda_nei*L
% M = cross-covariance matrix of size (3 x Nfeat)x(Nvox)
% Model weights are then selected based on optimal reg. param. pairs 
% See 'Pseudocode for SPIN-VM' from the paper for more details
X = single(X);
M = -X'*Y;
Mtilde = M*U;
Art = X'*X;
IA = eye(length(Art));
ureg_lap = unique(reg_lap);
ureg_reg = unique(reg_reg);
wsall = single(zeros(size(M,1),size(U,1),length(ureg_reg),length(ureg_lap)));
ws = single(zeros(size(X,2),size(U,1)));
for i=1:length(ureg_reg) % lambda_feat
    A=(Art+single(ureg_reg(i)*IA));% X'*X+lambda*I
    [Q, D] = eig(A); % Find eigenvalues and eigenvectors of X'*X+lambda*I
    D = diag(D);
    Qt = Q';
    sizA = size(A,1);
    Drep = repmat(D,1,size(U,1));
    for j=1:length(ureg_lap) % lambda_nei
        P = spinvm_rep(ureg_lap(j)*S,Drep,sizA);
        Xtilde = -Q*(P.*(Qt*Mtilde));
        w = Xtilde*U';
        w = real(w);
        w = single(w); % Save as single to reduce size
        wsall(:,:,i,j) = w; % Calculated model weights for each pair of reg. parameters
    end
end

for vx=1:size(U,1)
        [~,reg_ind,~] = intersect(ureg_reg,reg_reg(vx));
        [~,lap_ind,~] = intersect(ureg_lap,reg_lap(vx));
        ws(:,vx) = wsall(:,vx,reg_ind,lap_ind);
end
