function f = csindex(clust,X,DXX,varargin)
% CSINDEX Evaluation based on the CS criterion.
%   CSINDEX(CLUST, X, DXX) computes the CS criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the lowest index value.
%
%   CLUST is a numeric vector that represents a clustering solution.
%   X is an N-by-P data matrix with one row per observation and one
%   column per variable. DXX is an N-by-N dissimilarity matrix.
%
%   V = CSINDEX(CLUST, X) returns a positive numeric value corresponding to
%   the CS index. The dissimilarity matrix DXX is obtained
%   from X using the Euclidean distance.
%
%   V = CSINDEX(CLUST, X, DXX) returns a positive numeric value corresponding to
%   the CS index.
%
%   V = CSINDEX(..., 'DISTANCE', value) computes the CS index using
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
%   eva   = csindex(clust,meas,DXX);
%
%   See also EVALCVI, CVICONFIG, CVNNINDEX, COPINDEX, DUNNINDEX, DBINDEX
%
%   Reference:
%   ----------
%   Chien Hsing Chou, M. C. Su, and E. Lai, 
%   "A New Cluster Validity Measure and its Application to Image Compression," 
%   Pattern Analysis and Applications, 
%   Vol. 7, No. 2, pp. 205–220, 2004.
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
% Validation of the clustering solution 
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
num = 0;
for i = 1:K
   XC = DXX(clust==i,clust==i);
   num = num+mean(max(XC,[],2));
   
   members = (clust == cnames(i));
   M(i,:) = mean(X(members,:),1);
end

% Inter-cluster dispersion
Mij = feval(pfun,M',M');
Mij(logical(eye(K))) = Inf;

% CVI Computation
f = num/sum(min(Mij,[],2));