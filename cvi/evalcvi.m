function struct = evalcvi(clust, cvi, data, varargin)
% EVALCVI Evaluation of clustering solutions using a cluster validity index.
%   EVAL = EVALCVI(CLUST, CVI, DATA) creates an evaluation object
%   which can be used for estimating the number of clusters on data.
%
%   CLUST represents a set of clustering solutions. 
%   It must have N rows and contain integers. Column J contains the
%   cluster indices for each of the N points in the Jth clustering solution.
%   CVI is a string representing the criterion to be used. 
%   See the list of available CVIs by typing 'help cviconfig'.
%   DATA can be an N-by-P data matrix with one row per observation and one
%   column per variable (X) or it can be an N-by-N dissimilarity matrix (DXX).
%
%   If DATA is a feature set (X), it should be used when the CVI is:
%   'ch','db','xb','gd41','gd51','gd33','gd43','gd53','sdbw','pbm', 'cs',
%   'db2','sf','sym','sdb','sdi','cop','sv','wb','dbcv','lccv', or 'ssdd'.
%
%   If DATA is a dissimilarity matrix (DXX), it should be used when the CVI is:
%   'dunn','ci','sil','gd31','cvnn', or 'cvdd'
%
%
%   EVAL return an object with the following properties:
%   	OptimalK         - The optimal number of clusters suggested.
%       FitnessK         - The fitness value of the best clustering solution.
%       InspectedK       - List of the number of clusters inspected.
%       FitnessValues    - The fitness values for each number of clusters.
%
%   EVAL = EVALCVI(..., 'PARAM1', VALUE1, 'PARAM1', VALUE2) accepts one or two
%   comma-separated optional argument name/value pairs. Parameters are:
%
%   'Distance' -  perform the evaluation using a built-in measures:
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
%   Kmax = 6;
%   clust = zeros(size(meas,1),Kmax);
%   for k=1:Kmax
%       clust(:,k) = kmeans(meas,k,'distance','sqeuclidean');
%   end
%   eva = evalcvi(clust,'ch',meas);
%
%   See also CHINDEX, SILINDEX, DUNNINDEX
%
% ------------------------------------------------------------------------
%   Version 1.0 (Matlab R2020b Unix)
%   Copyright (c) 2021, A. Jose-Garcia and W. Gomez-Flores
% ------------------------------------------------------------------------

%Parameters validations
if nargin > 3
    [varargin{:}] = convertStringsToChars(varargin{:});
end

pnames = {'distance'}; pdvals = {'euc'};
[Dtype] = internal.stats.parseArgs(pnames, pdvals, varargin{:});

% Get the settings for the CVI
[cvifun,opt] = cviconfig(cvi);

% Number of clustering solutions in CLUST
nclust = size(clust,2);

Klist = NaN(1,nclust);
Values = NaN(1,nclust);

% Perform the evaluation of the cvi for each clustering solution
for i=1:nclust
    Klist(i) = numel(unique(clust(:,i)));
    Values(i) = feval(cvifun, clust(:,i), data);        
end

% Set NaN values if the cvi value is -inf or +inf
Values(Values==-inf | Values == inf) = NaN;

% Find the best number of clusters and fitness value
if strcmpi(opt,'mx')
    [bestFit,p] = max(Values);
else
    [bestFit,p] = min(Values);
end

% Return the evaluation object
struct.OptimalK = Klist(p);
struct.FitnessK = bestFit;
struct.InspectedK  = Klist;
struct.FitnessValues = Values;

end