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
addpath([pwd '/selection']);
addpath([pwd '/cvi']);
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
X = MODEData.X;
%% Run MODEC algorith
rng default; rng(1);
OUT = modec(MODEData);

%% Supervised Model Selection using ARI
[Tb,ARI_val,idx1,ARIvalues] = supervised(OUT.PClrs, MODEData.T);
disp(['Best ARI = ' num2str(ARI_val) ' | Best ID = ' num2str(idx1) ' | PFA size = ' num2str(size(OUT.PClrs,2))]);

[idx2] = unsupervised(OUT.PFront(1:end-1,:));
%% Figure 1: 2D-Scatter plot
figure(1)
dotp = 12; mg = [0.10 0.10 0.10];

pscat = scatter(X(:,1),X(:,2),dotp,'o','MarkerEdgeColor',mg,'MarkerFaceColor',mg);
pscat.MarkerFaceAlpha = .2;
pscat.MarkerEdgeAlpha = .5;

title('Unlabeled data');
niceplot(8);
set(gca,'xgrid','on','ygrid','on');

set(gcf, 'color','white');
set(gcf, 'renderer', 'painters');
figuresize(5,5,'centimeters');

%% Figure 2: Paret front plot
figure(2)
plotPFA(OUT.PFront,2,idx2)

title('Pareto front approximation');
xlabel('XB index'); ylabel('PBM index');
niceplot(8); set(gca,'xgrid','on','ygrid','on');

set(gcf, 'color','white');
set(gcf, 'renderer', 'painters');
figuresize(6,5,'centimeters');
%% Figure 3: Scatter plot
figure(3)
dotp = 12;
map = colormap(lines);

kT = unique(Tb);
for i = 1:numel(kT)
    idx = Tb==kT(i);
    pscat = scatter(X(idx,1),X(idx,2),dotp,'o','MarkerEdgeColor',map(i,:),'MarkerFaceColor',map(i,:));
    pscat.MarkerFaceAlpha = .2;
    pscat.MarkerEdgeAlpha = .5;
    hold on;
end
title('Data clustering');
txt = ['ARI = ' num2str(ARI_val,'%.2f')]; text(0.5,0.05,txt);

niceplot(8);
set(gca,'xgrid','on','ygrid','on');
%set(gca,'yTick',0:0.5:1,'xLim', [0 1]);

set(gcf, 'color','white');
set(gcf, 'renderer', 'painters');
figuresize(5,5,'centimeters');