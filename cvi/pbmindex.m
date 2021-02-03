function f = pbmindex(X,clrs,K,pfun)
global Xmean;

% Validacion del agrupamiento
N = numel(clrs);
clusts = unique(clrs);
Nk = accumarray(clrs,ones(N,1),[K,1]);

if numel(clusts) ~= K || sum(Nk<2)
    f = -inf;
    return;
end

% Dispersion intra-cluster
M = NaN(K,size(X,2));
sumD = zeros(K,1);
for i = 1:K
    members = (clrs == clusts(i));
    M(i,:) = mean(X(members,:),1);
    sumD(i)= sum(feval(pfun,X(members,:)',M(i,:)').^2);          
end
E1 = sum(feval(pfun,X',Xmean));
EK = sum(sumD);
E = E1/EK;

% Dispersion inter-cluster
Dij = feval(pfun,M',M');
D = max(Dij(~eye(K)));

% Evaluacion del IVG
K1 = 1/K;
f = (K1*E*D)^2;