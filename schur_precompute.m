function [U, S] = schur_precompute(sub,roiname,cube_size,mask_size)
% Precompute Schur decomposition of the graph Laplacian and save for further use

load(sprintf('./data/_roi_neighbor_%s_%s_%dc%d.mat',sub,roiname,cube_size,mask_size),'L')
[U,S] = schur(L,'complex');
S = diag(S);
sc = sprintf('./data/%sschurvars%s%dc%d.mat', sub,roiname,cube_size,mask_size);
eval(['save -v7.3 ' sc(1:end-4) '.mat U S ']);