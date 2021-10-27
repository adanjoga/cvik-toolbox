function f = symdunnindex(clust,X,DXX,varargin)
% SYMDBINDEX Evaluation based on the  Dunn based on Symmetry (SDunn) criterion.
%   SYMDUNNINDEX(CLUST, X, DXX) computes the Dunn based on Symmetry (SDunn) criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the highest index value.
%
%   X is an N-by-P data matrix with one row per observation and one
%   column per variable. CLUST is a numeric vector that represents 
%   a clustering solution. By default, the SDunn index uses the Symetry
%   distance and internally uses the Euclidean distance.
%
%   V = SYMDUNNINDEX(CLUST, X) returns a positive numeric value corresponding to
%   the SDunn index. The dissimilarity matrix DXX is obtained
%   from X using the Euclidean distance.
%
%   V = SYMDUNNINDEX(CLUST, X, DXX) returns a positive numeric value corresponding to
%   the SDunn index.
%
%   V = SYMDUNNINDEX(..., 'DISTANCE', value) computes the SDunn index using
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
%   eva   = symdunnindex(clust,meas,DXX);
%
%   See also EVALCVI, CVICONFIG, SILINDEX, DUNNINDEX
%
%   Reference:
%   ----------
%   S. Saha and S. Bandyopadhyay, "Performance Evaluation of Some 
%   Symmetry-Based Cluster Validity Indexes," 
%   IEEE Transactions on Systems, Man, and Cybernetics, Part C (Applications and Reviews), 
%   Vol. 39, No. 4, pp. 420–425, 2009.
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
    f = -inf;
    return;
end
% Evaluacion del IVG
M = NaN(K,size(X,2));
dmax = zeros(1,K);
dmin = zeros(K);
for i = 1:K
    members = (clust == cnames(i));
    M(i,:) = mean(X(members,:),1);
    
    dps = symdist(X(clust == i,:),M(i,:),Nk(i),pfun);
    dmax(i) = max(dps);    
    for j = setdiff(1:K,i)
        dmin(i,j) = min(min(DXX(clust == i,clust==j)));
    end
end
f = min(dmin(~eye(K))/max(dmax));