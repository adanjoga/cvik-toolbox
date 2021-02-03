function f = xbindex(X,clrs,K,pfun)
% Validacion del agrupamiento
N = numel(clrs);
clusts = unique(clrs);
Nk = accumarray(clrs,ones(N,1),[K,1]);

if numel(clusts) ~= K || sum(Nk<2)
    f = inf;
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
num = sum(sumD);

% Dispersion inter-cluster
den = feval(pfun,M',M').^2;
den(logical(eye(K))) = Inf;
den = min(den(:));

% Evaluacion del IVG
f = (1/N)*(num/den);

