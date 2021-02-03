function f = svindex(X,clrs,K,pfun)
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
maxD = zeros(K,1);
for i = 1:K
    members = (clrs == clusts(i));
    M(i,:) = mean(X(members,:),1);
    maxD(i) = max(feval(pfun,X(members,:)',M(i,:)'));
end
den = sum(maxD,1);

% Dispersion inter-cluster
Cij = feval(pfun,M',M');
Cij(logical(eye(K))) = Inf;
num = sum(min(Cij,[],1));

f = num/den;

