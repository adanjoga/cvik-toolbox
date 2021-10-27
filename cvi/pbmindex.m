function f = pbmindex(clust,X,varargin)
% PBMINDEX Evaluation based on the PBM criterion.
%   SFINDEX(CLUST, X) computes the PBM criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the highest index value.
%
%   X is an N-by-P data matrix with one row per observation and one
%   column per variable. CLUST is a numeric vector that represents 
%   a clustering solution. By default, the PBM index uses
%   the Euclidean distance between points in X.
%
%   V = PBMINDEX(CLUST, X) returns a positive numeric value corresponding to
%   the PBM index.
%   
%   V = PBMINDEX(..., 'DISTANCE', value) computes the PBM index using
%   a specified distance measure. The available built-in measures are:
%       'euc'           - Euclidean distance (the default).
%       'neuc'          - Normalized Euclidean distance.
%       'cos'           - Cosine similarity.
%       'pcorr'         - Pearson's correlation coefficient.
%       'scorr'         - Spearman's correlation coefficient.
%       'lap'           - Laplacian distance.
%
%
%   Example:
%   -------
%   load fisheriris;
%   clust = kmeans(meas,3,'distance','sqeuclidean');
%   eva   = pbmindex(clust, meas);
%
%
%   See also EVALCVI, CVICONFIG, CHINDEX, XBINDEX, DBINDEX, SFINDEX, DBCVINDEX
%
%
%   Reference:
%   ----------
%   U.  Maulik  and  S.  Bandyopadhyay, "Performance  Evaluation  of  Some  
%   Clustering Algorithms  and  Validity  Indices," 
%   IEEE  Transactions  on  Pattern  Analysis  and  Machine  Intelligence, 
%   Vol. 24, No. 12, pp. 1650–1654, 2002.
%
% ------------------------------------------------------------------------
%   Version 1.0 (Matlab R2020b Unix)
%   Copyright (c) 2021, A. Jose-Garcia and W. Gomez-Flores
% ------------------------------------------------------------------------

% Parameter validation
if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end
pnames = {'distance'}; pdvals = {'euc'};
[Dtype] = internal.stats.parseArgs(pnames, pdvals, varargin{:});
pfun = proxconfig(Dtype);
% ------------------------------------------------------------------------
% WARNING: If possible compute the Xmean value before calling the SF index
% global Xmean;
Xmean = mean(X,1)';  % Global mean value of input data X
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

% Intra-cluster cohesion (compactness)
M = NaN(K,size(X,2));
sumD = zeros(K,1);
for i = 1:K
    members = (clust == cnames(i));
    M(i,:) = mean(X(members,:),1);
    sumD(i)= sum(feval(pfun,X(members,:)',M(i,:)').^2);          
end
E1 = sum(feval(pfun,X',Xmean));
EK = sum(sumD);
E = E1/EK;

% Dispersion inter-cluster
Dij = feval(pfun,M',M');
D = max(Dij(~eye(K)));

% CVI Computation
K1 = 1/K;
f = (K1*E*D)^2;