% -------------------------------------------------------------------------
%% Multiple Distances Project
% Name:     Adan Jose Garcia
% Date:     Feb 18th, 2019

% The K-medoid algorithm
% -------------------------------------------------------------------------
%% Load variables and directories
clearvars; close all; clc;

addpath([pwd '/proximity']);
addpath([pwd '/datasets/synthetic']); addpath([pwd '/datasets/MED_matrices']);
addpath([pwd '/clustering']);
addpath([pwd '/cvi']);
addpath([pwd '/validation']);
addpath([pwd '/utils']);
%% variables regarding the clustering problem
iDS = 4;
DSnames = { 'Orange','Data_4_3','Data_6_2','R15','Twenty',...
            'TwoDiamonds','Square1','Sizes5','Data_5_2','Data_9_2',...
            'Part2','Inside','Spirals2','Chainlink','Atom',...
            'Flame','Ringauss','Spiralsquare','Spiralsizes5','Multidist'};
        
DSnames = {'Spiralsdata52','Spiralsdata92','Spiralsflame','Flamesize5'};
        
Params.dist = 'med_euc';    % {euclidean, med_euc , cosine}

% DMat = load('DMat_4dist.mat'); % load EXX, GEXX, CXX and GCXX
% nDS  = numel(DSnames);
%% Perform the experiment [Dataset]
disp(['Computing dataset: ' DSnames{iDS}]);
Data    = load(DSnames{iDS});
X       = Data.data(:,1:end-1);     % Dataset
T       = Data.data(:,end);         % True labeling of the dataset
Ktrue   = numel(unique(T));         % True number of clusters

if strcmpi(Params.dist,'euclidean')
    %DXX = DMat.EXX{iDS};
elseif strcmpi(Params.dist,'cosine')
    %DXX = DMat.CXX{iDS};
elseif strcmpi(Params.dist,'med_euc')
    DXmed = load([DSnames{iDS} '_euc']); % MED distance based on Euclidean
    DXX = DXmed.DXX;     
end

%% Run PAM algorithm
if strcmpi(Params.dist,'med_euc')
    cellout = kmedoids4(DXX,Ktrue);
    Yb = cellout{3};
    wcd = cellout{1};
    Med = cellout{4};
    M = X(Med,:)

else
    [Yb,M,wcd] = kmeans(X,Ktrue,'empty','singleton','replicate',5,'Distance','sqeuclidean');
end

ARI = pairwiseindex(T,Yb)

plotclusters2(X,T,Yb,M);
%% Decodification of solutions
% meds    = clr2med(T,Ktrue,DXX);
% Clr     = med2clr(meds.medoids,DXX);
% 
% ARIb = pairwiseindex(T,Clr)