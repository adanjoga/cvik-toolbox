function f = sdbwindex(X,clrs,K,pfun)
global sXX;

% Validacion del agrupamiento
N = numel(clrs);
clusts = unique(clrs);
Nk = accumarray(clrs,ones(N,1),[K,1]);
if numel(clusts) ~= K || sum(Nk<2)
    f = inf;
    return;
end

% Intra-cluster variance
sM  = zeros(K,size(X,2)); % Varianza del cluster Ci
for i = 1:K
    sM(i,:) = var(X(clrs==i,:),0,1);
end
sMM = sqrt(sum(sM.^2,2)); % 2-norma |x| = (x.x)^0.5
% Varianza del dataset: sX  = var(X,0,1);
% 2-norma del dataset:  sXX = sqrt(sum(sX.^2,2));
Scat = (1/K)*sum(sMM/sXX);

% Average standard deviation of clusters
stdev = (1/K)*sqrt(sum(sMM));

% Inter-cluster density
M = NaN(K,size(X,2));
Dmin = zeros(N,1);
for i = 1:K
    members = (clrs == clusts(i));
    M(i,:) = mean(X(members,:),1);
    Dmin(members) = feval(pfun,X(members,:)',M(i,:)');
end

aux1 = zeros(1,K);
for i = 1:K
   c = 0;
   Di = density(Dmin(clrs==i,:),[],stdev,pfun);
   aux2 = zeros(1,K-1);
   for j = setdiff(1:K,i)
       c = c+1;
       Uij = 0.5*(M(i,:)+M(j,:));
       Dij = density(X(clrs==i|clrs==j,:),Uij,stdev,pfun);
       Dj  = density(Dmin(clrs==j,:),[],stdev,pfun);
       aux2(c) = Dij/max(max(Di,Dj),1);
   end
   %aux2(isinf(aux2)|isnan(aux2)) = 0;
   aux1(i) = sum(aux2);
end
Dens_bw = (1/(K*(K-1)))*sum(aux1);

% Evaluacion del IVG
f = Scat + Dens_bw;

if isnan(f)
    inf;
end
%***********************************************************************
function D = density(Z,U,sd,pfun)
if isempty(U)
    D = sum(Z <= sd,1);
else
    D = sum(feval(pfun,Z',U') <= sd,1);
end