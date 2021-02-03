function f = cindex(X,clrs,K,pfun)
global DXX;
% Validacion del agrupamiento
N = numel(clrs);
clusts = unique(clrs);
Nk = accumarray(clrs,ones(N,1),[K,1]);
if numel(clusts) ~= K || sum(Nk<2)
    f = Inf;
    return;
end
% DXX = feval(pfun,X',X');
idG = triu(~eye(N),1);
A = DXX(idG);
As = sort(A);
S = 0;
for i = 1:K
   id1 = double(clrs==i);
   id2 = logical(id1*id1');
   B = id2(idG);
   S = S+sum(A(B));
end
nw = (sum(Nk.^2)-N)/2;
nt = N*(N-1)/2;
Smin = sum(As(1:nw));
Smax = sum(As(nt-nw+1:nt));
f = (S-Smin)/(Smax-Smin);