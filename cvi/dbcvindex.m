function valid = dbcvindex(clust,X,varargin)
% DBCVINDEX Evaluation based on the Density-Based criterion.
%   DBCVINDEX(CLUST, X) computes the Density-Based criterion 
%   value which can be used for estimating the number of clusters on data.
%   The optimal number of clusters is the solution with the highest index value.
%
%   X is an N-by-P data matrix with one row per observation and one
%   column per variable. CLUST is a numeric vector that represents 
%   a clustering solution. By default, the Density-Based index uses
%   the Euclidean distance between points in X.
%
%   V = DBCVINDEX(CLUST, X) returns a positive numeric value corresponding to
%   the Density-Based index.
%   
%   V = DBCVINDEX(..., 'DISTANCE', value) computes the Davies-Bouldin index using
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
%   eva   = dbcvindex(clust, meas);
%
%
%   See also EVALCVI, CVICONFIG, CHINDEX, XBINDEX, PBMINDEX, SFINDEX, DBINDEX
%
%
%   Reference:
%   ----------
%   D. Moulavi, P. A. Jaskowiak, R. J. G. B. Campello, A. Zimek, and J. Sander.
%   "Density-Based Clustering Validation".  
%   In SIAM International Conference on Data Mining, 
%   pp. 839–847, Philadelphia, PA, 2014.
%   NOTE: Code provided by the authors
% ------------------------------------------------------------------------
%   Version 1.0 (Matlab R2020b Unix)
%   Copyright (c) 2021, A. Jose-Garcia and W. Gomez-Flores
% ------------------------------------------------------------------------
%Parameter validation
if nargin > 2
    [varargin{:}] = convertStringsToChars(varargin{:});
end
pnames = {'distance'}; pdvals = {'euc'};
[mydist] = internal.stats.parseArgs(pnames, pdvals, varargin{:});
%pfun = proxconfig(mydist);

% ------------------------------------------------------------------------
% Validation of the clustering solution 
clusters  = unique(clust);
DXX      = squareform(pdist(X)).^2;
%treating singleton clusters
for i=1:length(clusters)
    if(sum(clust == clusters(i)) == 1)
        clust(clust == clusters(i)) = 0;
        clusters(i) = 0;
    end
end

%clusters without noise and singletons
clusters  = setdiff(clusters,0);

if (isempty(clusters) || (length(clusters) == 1))
    valid = 0;
    return;
end

X      = X(clust~=0,:);
DXX      = DXX(clust~=0,clust~=0);

poriginal = clust;
clust = clust(clust~=0);


nclusters = length(clusters);

[nobjects nfeatures] = size(X);

d_ucore_cl = zeros(1,nobjects);
compcl     = zeros(1,nclusters);
int_edges  = cell(1,nclusters);
int_node_data   = cell(1,nclusters);
for i=1:nclusters
    
    objcl  = find(clust == clusters(i));
    
    nuobjcl = length(objcl);
            
    [d_ucore_cl(objcl) mr] = matrix_mutual_reachability_distance(nuobjcl, DXX(objcl,objcl),nfeatures);  % ucore distance of each object in its own cluster
     
    G.no_vertices = nuobjcl;
    G.MST_edges   = zeros(nuobjcl-1,3);
    G.MST_degrees = zeros(nuobjcl,1);
    G.MST_parent  = zeros(nuobjcl,1);
    
    [Edges Degrees] = MST_Edges(G, 1,mr);
    
    int_node     = find(Degrees~=1);
    int_edg1     = find(ismember(Edges(:,1),int_node));    
    int_edg2     = find(ismember(Edges(:,2),int_node));    
    int_edges{i} = intersect(int_edg1,int_edg2);
    
    
    if (~isempty(int_edges{i}))
        compcl(i) = max(Edges(int_edges{i},3));
    else 
        compcl(i) = max(Edges(:,3));         
    end
    int_node_data{i} = objcl(int_node);
    if isempty(int_node_data{i})
        int_node_data{i} = objcl;
    end
end

sep_point = zeros(nobjects,nobjects);
for i=1:(nobjects-1)
    for j=i:nobjects
        sep_point(i,j) = max([DXX(i,j) d_ucore_cl(i) d_ucore_cl(j)]);        
        sep_point(j,i) = sep_point(i,j);
    end
end

valid = 0;
sepcl = Inf(nclusters,1);
for i=1:nclusters
    other_cls = setdiff(clusters,clusters(i)); 
    
    sep = zeros(1,length(other_cls));
    for j=1:length(other_cls)                
        sep(j) = min(min(sep_point(int_node_data{i},int_node_data{clusters==other_cls(j)}))); 
    end
    sepcl(i) = min(sep);
    dbcvcl = (sepcl(i) - compcl(i)) / max(compcl(i),sepcl(i));
    valid = valid + (dbcvcl * sum(clust == clusters(i)));
end

valid = valid / length(poriginal); 

end    
