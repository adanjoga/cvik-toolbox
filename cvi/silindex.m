function f = silindex(clust, DXX)
% SILINDEX Evaluation based on the Silhouette index.
%   SILINDEX(CLUST, DXX) computes the Silhouette criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the highest index value.
%
%   DXX is an N-by-N dissimilarity matrix.
%   CLUST is a numeric vector that represents a clustering partition. 
%
%   V = SILINDEX(CLUST, DXX) returns a positive numeric value corresponding to
%   the Silhouette index.
%   
%
%   Example:
%   -------
%   load fisheriris;
%   clust = kmeans(meas,3,'distance','sqeuclidean');
%   DXX = pdist2(meas,meas,'Euclidean');
%   eva   = silindex(clust, DXX);
%
%   See also EVALCVI, CVICONFIG, CVNNINDEX, CVDDINDEX, DUNNINDEX, CINDEX
%
%
%   Reference:
%   ----------
%   Peter J. Rousseeuw, "Silhouettes:  A Graphical Aid to the Interpretation 
%   and Validation of Cluster Analysis," 
%   Journal of Computational and Applied Mathematics, 
%   Vol. 20, pp. 53-65, 1987.
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
DX = DXX.^2;
S = zeros(1,K);
for i = 1:K
    XC = DX(clust==cnames(i),clust==cnames(i));
    ai = sum(XC,2)/max(Nk(i)-1, 1); % modification 25-10-16
    %ai = mean(XC,2);
    dd = zeros(Nk(i),K-1);
    k = 0;
    for j = setdiff(1:K,i)
        k = k+1;
        XK = DX(clust==cnames(i),clust==cnames(j));
        dd(:,k) = mean(XK,2);
    end
    bi = min(dd,[],2);
    S(i) = sum((bi-ai)./max([ai bi],[],2),1);
end
f = (1/N)*sum(S);
