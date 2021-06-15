function f = dbindex(X,clrs,distance)
% Validacion del agrupamiento
N = numel(clrs);
clusts = unique(clrs);
K = numel(clusts);
Nk = accumarray(clrs,ones(N,1),[K,1]);

if numel(clusts) ~= K || sum(Nk<2)
    f = inf;
    return;
end

if nargin < 3 || isempty(distance)
    distType = 'euc';
else
    distType = distance;
end
pfun = proxconfig(distType);

% ------------------------------------
% Dispersion intra-cluster
M = NaN(K,size(X,2));
sumD = zeros(K,1);
for i = 1:K
    members = (clrs == clusts(i));
    M(i,:) = mean(X(members,:),1);
    sumD(i)= sum(feval(pfun,X(members,:)',M(i,:)'));          
end
Si = sumD./Nk;

% Dispersion inter-cluster
Mij = feval(pfun,M',M');
R = -inf(K);
for i = 1:K
    for j = setdiff(1:K,i)
        R(i,j) = (Si(i)+Si(j))/Mij(i,j);
    end
end
% Evaluacion del IVG
f = mean(max(R,[],2));
