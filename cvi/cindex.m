function f = cindex(clust, DXX)
% CINDEX Evaluation based on the C-index.
%   CINDEX(CLUST, DXX) computes the C-criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the lowest index value.
%
%   DXX is an N-by-N dissimilarity matrix.
%   CLUST is a numeric vector that represents a clustering partition. 
%
%   V = CINDEX(CLUST, DXX) returns a positive numeric value corresponding to
%   the C-index.
%   
%   Example:
%   -------
%   load fisheriris;
%   clust = kmeans(meas,3,'distance','sqeuclidean');
%   DXX = pdist2(meas,meas,'Euclidean');
%   eva   = cindex(clust, DXX);
%
%   See also EVALCVI, CVICONFIG, SILINDEX, CVNNINDEX, CVDDINDEX, DUNNINDEX
%
%
%   Reference:
%   ----------
%   L. J. Hubert and J. R. Levin, "A General Statistical Framework for 
%   Assessing Categorical Clustering in Free Recall," 
%   Psychological Bulleting, 
%   Vol. 83, pp. 1072–1080, 1976.
%
% ------------------------------------------------------------------------
%   Version 1.0 (Matlab R2020b Unix)
%   Copyright (c) 2021, A. Jose-Garcia and W. Gomez-Flores
% ------------------------------------------------------------------------

% Validation of the clustering solution 
N = numel(clust);
cnames = unique(clust);
K = numel(cnames);
Nk = accumarray(clust,ones(N,1),[K,1]);
if sum(Nk<1) || K==1
    f = Inf;
    return;
end

% CVI Computation
idG = triu(~eye(N),1);
A = DXX(idG);
As = sort(A);
S = 0;
for i = 1:K
   id1 = double(clust==i);
   id2 = logical(id1*id1');
   B = id2(idG);
   S = S+sum(A(B));
end
nw = (sum(Nk.^2)-N)/2;
nt = N*(N-1)/2;
Smin = sum(As(1:nw));
Smax = sum(As(nt-nw+1:nt));
f = (S-Smin)/(Smax-Smin);