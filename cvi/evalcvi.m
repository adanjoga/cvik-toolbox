function struct = evalcvi(X, clust, cvi, distance)
% EVALCVI Evaluation of clustering solutions using a cluster validity index.
%   EVAL = EVALCVI(X, CLUST, CVI) creates an evaluation object
%   which can be used for estimating the number of clusters on data.
%
%   X must be an N-by-P matrix of data with one row per observation and one
%   column per variable. CLUST is a numeric matrix representing a set of  
%   clustering solutions. It must have N rows and contain integers. 
%   Column J contains the cluster indices for each of the N points for the 
%   J-th clustering solution. CVI is a string representing the criterion
%   to be used. See the list of available CVIs by typing 'help cviconfig'.
%
%   EVAL = EVALCVI(X, CLUST, CVI, DISTANCE) perform the evaluation using
%   the specified distance measure. The available built-in measures are:
%       'euc'           - Euclidean distance (the default).
%       'neuc'          - Normalized Euclidean distance.
%       'cos'           - Cosine similarity.
%       'pcorr'         - Pearson's correlation coefficient.
%       'scorr'         - Spearman's correlation coefficient.
%       'lap'           - Laplacian distance.
%
%   EVALCVI return an object with the following properties:
%   	OptimalK         - The optimal number of clusters suggested.
%       FitnessK         - The fitness value of the best clustering solution.
%       InspectedK       - List of the number of clusters inspected.
%       FitnessValues    - The fitness values for each number of clusters.
%
%
%   Example:
%   -------
%   load fisheriris;
%   Kmax = 6;
%   clust = zeros(size(meas,1),Kmax);
%   for k=1:Kmax
%       clust(:,k) = kmeans(meas,k,'distance','sqeuclidean');
%   end
%   eva = evalcvi(meas,clust,'ch');
%
%
%   See also CHINDEX, SILINDEX, DUNNINDEX
%
% ------------------------------------------------------------------------
%   Version 1.0 (Matlab R2020b Unix)
%   Copyright (c) 2021, A. Jose-Garcia and W. Gomez-Flores
% ------------------------------------------------------------------------

if nargin < 4 || isempty(distance)
    distType = 'euc';
else
    distType = distance;
end

% Get the settings for the CVI
[cvifun,opt] = cviconfig(cvi);

% Number of clustering solutions in CLUST
nclust = size(clust,2);

Klist = NaN(1,nclust);
Values = NaN(1,nclust);

% Perform the evaluation of the cvi for each clustering solution
for i=1:nclust
    Klist(i) = numel(unique(clust(:,i)));
    Values(i) = feval(cvifun, X, clust(:,i), distType);        
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