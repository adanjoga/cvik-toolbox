function f = copindex(clust,X,DXX,varargin)
% COPINDEX Evaluation based on the COP criterion.
%   COPINDEX(CLUST, X, DXX) computes the COP criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the lowest index value.
%
%   CLUST is a numeric vector that represents a clustering solution.
%   X is an N-by-P data matrix with one row per observation and one
%   column per variable. DXX is an N-by-N dissimilarity matrix.
%
%   V = COPINDEX(CLUST, X) returns a positive numeric value corresponding to
%   the SDunn index. The dissimilarity matrix DXX is obtained
%   from X using the Euclidean distance.
%
%   V = COPINDEX(CLUST, X, DXX) returns a positive numeric value corresponding to
%   the COP index.
%
%   V = COPINDEX(..., 'DISTANCE', value) computes the COP index using
%   a specified distance measure. The available built-in measures are:
%       'euc'           - Euclidean distance (the default).
%       'neuc'          - Normalized Euclidean distance.
%       'cos'           - Cosine similarity.
%       'pcorr'         - Pearson's correlation coefficient.
%       'scorr'         - Spearman's correlation coefficient.
%       'lap'           - Laplacian distance.
%
%   Example:
%   -------
%   load fisheriris;
%   clust = kmeans(meas,3,'distance','sqeuclidean');
%   DXX = pdist2(meas,meas,'Euclidean');
%   eva   = copindex(clust,meas,DXX);
%
%   See also EVALCVI, CVICONFIG, CVNNINDEX, CSINDEX, DUNNINDEX, DBINDEX
%
%   Reference:
%   ----------
%   Ibai Gurrutxaga, et al., "SEP/COP: An Efficient Method to Find the 
%   Best Partition in Hierarchical Clustering Based on a New Cluster Validity Index," 
%   Pattern Recognition, 
%   Vol. 43, No. 10, pp. 3364–3373, 2010.
%
% ------------------------------------------------------------------------
%   Version 1.0 (Matlab R2020b Unix)
%   Copyright (c) 2021, A. Jose-Garcia and W. Gomez-Flores
% ------------------------------------------------------------------------
%Parameter validation
if nargin > 3
    [varargin{:}] = convertStringsToChars(varargin{:});
end
pnames = {'distance'}; pdvals = {'euc'};
[Dtype] = internal.stats.parseArgs(pnames, pdvals, varargin{:});
pfun = proxconfig(Dtype);

if nargin < 3
    DXX = real(feval(pfun,X',X'));
end
% ------------------------------------------------------------------------
%Validation of the clustering solution 
N = numel(clust);
cnames = unique(clust);
K = numel(cnames);
Nk = accumarray(clust,ones(N,1),[K,1]);
if sum(Nk<1) || K==1
    f = inf;
    return;
end

% Intra-cluster cohesion (compactness)
M = NaN(K,size(X,2));
sumD = zeros(K,1);
for i = 1:K
    members = (clust == cnames(i));
    M(i,:) = mean(X(members,:),1);
    sumD(i)= sum(feval(pfun,X(members,:)',M(i,:)').^2);          
end
intra = (sumD./Nk)';

% Inter-cluster dispersion
inter = zeros(1,K);
for i = 1:K
    dd = zeros(Nk(i),K-1);
    k = 0;
    for j = setdiff(1:K,i)
        k = k+1;
        XC = DXX(clust==i,clust==j);
        dd(:,k) = max(XC,[],2);
    end
    inter(i) = min(dd(:));
end
% CVI Computation
f = (1/N)*sum(Nk.*(intra./inter)');
