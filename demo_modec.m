%% DEMO: MODEC algorithm appplied to the automatic clustering problem
% ------------------------------------------------------------------------
% The Cluster Validity Index Tooolbox (CVIT) for automatic determination 
% of clusters from clustering solution contains more than 70 functions (m-files)
% This toolbox was developed with MATLAB R2014a.
%
% Two CVIs are used as function objectives
%
% IMPORTANT: First run "RUN_ME_FIRST.m" file to add this toolbox to search path.
%------------------------------------------------------------------------
clc; clear all; close all;

addpath([pwd '/clustering']);
addpath([pwd '/proximity']);
addpath([pwd '/datasets']);
addpath([pwd '/validation']);
addpath([pwd '/cvi']);
addpath([pwd '/selection']);
addpath([pwd '/utils']);

% List of available cluster validity indices (CVIs)
CVInames = {'xb','ch','sf','pbm','cs',...
            'gd31','gd41','gd51','gd33','gd43',...
            'gd53','db2','db','cop','sil',...
            'dunn','sv','sym','sdunn','sdb',...
            'sdbw','cind'};

% List of available distances
Distnames = {'euc','neuc','cos','pcorr','scorr','lap'};

% List of datasets provided
DSnames = {'Data_4_3','Data_5_2','Sizes5', 'Iris'}; 

% ------------------------------------------------------------------------
%% Variables regarding the optimization problem
dataID = 3;
% load fisheriris;
% X = meas; T = [ones(50,1); ones(50,1)*2; ones(50,1)*3]; % Iris labels, k=3
% MODEData.X = X; MODEData.T=T;

Data = load(DSnames{dataID});         % Load the dataset listed in DS_names[]
MODEData.X = Data.data(:,1:end-1);     % Dataset
MODEData.T = Data.data(:,end);         % True labeling of dataset

MODEData.KMAX = 10;                 % Maximum number of clusters
MODEData.NOBJ = 2;                  % Number of objectives          
MODEData.CVIs = {'XB','PBM'};    % XB_PBM, XB_DB

MODEData.POPSIZE = 50;          % Population size
MODEData.MAXGEN = 500;          % Generation bound 

%% Run MODEC algorith
%rng default; rng(0);
OUT = modec(MODEData);

%% Supervised Model Selection using ARI
[Clr,ARIb,idx1,ARIvalues] = supervised(OUT.PClrs, MODEData.T);
disp(['Best ARI = ' num2str(ARIb) ' | Best ID = ' num2str(idx1) ' | PFA size = ' num2str(size(OUT.PClrs,2))]);

plotPFA(OUT.PFront,2,idx1)