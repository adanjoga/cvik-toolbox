%% DEMO: Evaluation of clustering solutions generated by k-means 
% ------------------------------------------------------------------------
% The Cluster Validity Index Tooolbox (CVIT) for automatic determination 
% of clusters from clustering solution contains more than 70 functions (m-files)
% This toolbox was developed with MATLAB R2014a.
%
% IMPORTANT: First run "RUN_ME_FIRST.m" file to add this toolbox to search path.
%------------------------------------------------------------------------
clc; clear all; close all;

addpath([pwd '/proximity']);
addpath([pwd '/cvi']);
addpath([pwd '/datasets']);
addpath([pwd '/validation']);
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
DSnames = {'Data_4_3','Data_5_2','Moon', 'Iris', 'Spirals'};

% ------------------------------------------------------------------------
%% Evaluate a set of clustering solutions generated by k-means algorithm
cvi = CVInames{15};

OriData = load(DSnames{5}); Data = load('Spirals_med_euc.mat');
DXX = Data.DXX; % MED-based dissimilarity matrix
X = OriData.data(:,1:end-1); T = OriData.data(:,end);

% Generation of clusterings using k-means and varying the number of clusters
rng default; rng(0);
Kmax = 10;
clust = zeros(size(DXX,1),Kmax);
for k=1:Kmax
     clust(:,k) = kmedoids(DXX,k);
end

% Evaluation of the clustering solutions using the 'sil' index
eva = evalcvi(DXX, clust, cvi, 'DataType','relational');

Kb = eva.OptimalK; % Estimated k value
Tb = clust(:,Kb);  % Estimated clustering solution
ARI_val = pairwiseindex(T,Tb);
% ------------------------------------------------------------------------
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

% ------------------------------------------------------------------------
%% Figure 2: Line plot (CVI convergence)
figure(2)    
mg = [0.60 0.60 0.60]; mb =[0 0.4470 0.7410];

k_vals   = 1:Kmax;
cvi_vals = eva.FitnessValues; cvi_vals(1) = 0;
plot(k_vals,cvi_vals,'-o','Color',mb,'MarkerSize',3,'LineWidth',1.5,'MarkerFaceColor',mb); hold on;

title([upper(cvi) ': Optimal number of clusters']);
xline(Kb,'Color',mg,'LineStyle','-.','LineWidth',1.5);
ylabel([upper(cvi) ' value']); xlabel('number of clusters');

niceplot(8); box on; 
set(gca,'YGrid', 'on'); set(gca,'xGrid', 'on'); 
set(gca,'xTick',1:1:Kmax, 'XLim',[1 Kmax]);
set(gcf, 'color','white');
set(gcf, 'renderer', 'painters');
figuresize(8,5,'centimeters');

% ------------------------------------------------------------------------
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
txt = ['ARI = ' num2str(ARI_val,'%.1f')]; text(0.4,0.1,txt);

niceplot(8);
set(gca,'xgrid','on','ygrid','on');
%set(gca,'yTick',0:0.5:1,'xLim', [0 1]);

set(gcf, 'color','white');
set(gcf, 'renderer', 'painters');
figuresize(5,5,'centimeters');