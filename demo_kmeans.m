%% DEMO: Evaluation of clustering solutions generated by k-means 
% ------------------------------------------------------------------------
% The Cluster Validity Index Tooolbox (CVIT) for automatic determination 
% of clusters from clustering solution contains more than 70 functions (m-files)
% This toolbox was developed with MATLAB R2014a.
%
% Developed by
%   Adan Jose-Garcia (adan.jose@cinvestav.mx)
%   Wilfrido Gomez Flores (wgomez@cinvestav.mx)
%
% IMPORTANT: First run "RUN_ME_FIRST.m" file to add this toolbox to search path.
%------------------------------------------------------------------------
clc; clear all; close all;

addpath([pwd '/proximity']);
addpath([pwd '/cvi']);
addpath([pwd '/datasets']);

% List of available cluster validity indices (CVIs)
CVInames = {'xb','ch','sf','pbm','cs',...
            'gd31','gd41','gd51','gd33','gd43',...
            'gd53','db2','db','cop','sil',...
            'dunn','sv','sym','sdunn','sdb',...
            'sdbw','cind'};

% List of available distances
Distnames = {'euc','neuc','cos','pcorr','scorr','lap'};
   
% List of some datasets provided
DSnames = {'Data_4_3','Data_5_2','Moon', 'Iris'};

% ------------------------------------------------------------------------
%% Evaluate a clustering solution generated by the k-means algorithm

% Example 1
load Data_4_3;
X = data(:,1:end-1);
clust = kmeans(X,4,'distance','sqeuclidean');
eva1 = chindex(clust,X);
%E2 = evalclusters(X,clust,'CalinskiHarabasz');

% Example 2
load fisheriris;
clust = kmeans(meas,3,'distance','sqeuclidean');
eva2 = chindex(clust,meas);

%% Evaluate a set of clustering solutions generated by k-means algorithm

% Load the Iris dataset
load fisheriris;

% Generation of clusterings using k-means and varying the number of clusters
Kmax = 6;
clust = zeros(size(meas,1),Kmax);
for k=1:Kmax
    clust(:,k) = kmeans(meas,k,'distance','sqeuclidean');
end

% Evaluation of the clustering solutions using the 'ch' index
eva = evalcvi(clust, 'ch', meas)

