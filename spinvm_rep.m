function [P] = spinvm_rep (S,Drep,sizA)
% This little script computes matrix P, of size (3 x Nfeat)x(Nvox), which
% is a key component in the computation of the final model weights
S = S';
Srep = repmat(S,sizA,1);
P = 1./(Drep+Srep);
end