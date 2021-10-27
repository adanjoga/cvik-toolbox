function f = gd31index(clust, DXX)
% GD31INDEX Evaluation based on the Generalized Dunn-31 index.
%   GD31INDEX(CLUST, DXX) computes the Generalized Dunn-31 criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the highest index value.
%
%   DXX is an N-by-N dissimilarity matrix.
%   CLUST is a numeric vector that represents a clustering partition. 
%
%   V = GD31INDEX(CLUST, DXX) returns a positive numeric value corresponding to
%   the Generalized Dunn-31 index.
%   
%
%   Example:
%   -------
%   load fisheriris;
%   clust = kmeans(meas,3,'distance','sqeuclidean');
%   DXX = pdist2(meas,meas,'Euclidean');
%   eva   = gd31index(clust, DXX);
%
%   See also EVALCVI, CVICONFIG, SILINDEX, CVNNINDEX, CVDDINDEX, DUNNINDEX
%
%   Reference:
%   ----------
%   James C. Bezdek and N. R. Pal, "Some New Indexes of Cluster Validity," 
%   IEEE Transactions on Systems, Man, and Cybernetics, Part B (Cybernetics), 
%   Vol. 28, No. 3, pp. 301–315, 1998.
%
% ------------------------------------------------------------------------
%   Version 1.0 (Matlab R2020b Unix)
%   Copyright (c) 2021, A. Jose-Garcia and W. Gomez-Flores
% ------------------------------------------------------------------------

N = numel(clust);
cnames = unique(clust);
K = numel(cnames);
Nk = accumarray(clust,ones(N,1),[K,1]);
if sum(Nk<1) || K==1
    f = -inf;
    return;
end

% CVI Computation
num = inf(K);
den = zeros(1,K);
for i = 1:K
    den(i) = max(max(DXX(clust==i,clust==i)));
    for j = setdiff(1:K,i)
        num(i,j) = mean2(DXX(clust==i,clust==j));
    end
end
f = min(num(:)/max(den));