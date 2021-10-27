function f = dunnindex(clust, DXX)
% DUNNINDEX Evaluation based on the Dunn index.
%   DUNNINDEX(CLUST,DXX) computes the Dunn criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the highest index value.
%
%   DXX is an N-by-N dissimilarity matrix.
%   CLUST is a numeric vector that represents a clustering partition. 
%
%   V = DUNNINDEX(CLUST,DXXX) returns a positive numeric value corresponding to
%   the Dunn index.   
%
%   Example:
%   -------
%   load fisheriris;
%   clust = kmeans(meas,3,'distance','sqeuclidean');
%   DXX = pdist2(meas,meas,'Euclidean');
%   eva   = dunnindex(clust,DXX);
%
%   See also EVALCVI, CVICONFIG, SILINDEX, CVNNINDEX, CVDDINDEX, CINDEX
%
%
%   Reference:
%   ----------
%   J. C. Dunn, "A fuzzy relative of the ISODATA process and its use in detecting 
%   compact well-separated clusters," 
%   Journal of Cybernetics, 
%   Vol. 3, No. 3, pp. 32–57, 1973.
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
    f = -inf;
    return;
end

% CVI Computation
dmin = inf(K);
dmax = zeros(1,K);
for i = 1:K
    dmax(i) = max(max(DXX(clust==i,clust==i)));
    for j = setdiff(1:K,i)
        dmin(i,j) = min(min(DXX(clust==i,clust==j)));
    end
end
f = min(dmin(:)/max(dmax));