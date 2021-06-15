function f = silindex(X,Xtype,clust,distance)

% Validation of the clustering solution 
N = numel(clust);
cnames = unique(clust);
K = numel(cnames);
Nk = accumarray(clust,ones(N,1),[K,1]);

if sum(Nk<1) || K==1
    f = -inf;
    return;
end

if nargin < 4 || isempty(distance)
    distance = 'euc';
end

if strcmp(Xtype,'feature')
    pfun    = proxconfig(distance);     % Defines the proximity measure function.
    DXX     = real(feval(pfun,X',X'));  % distance matrix based on pfun.
elseif strcmp(Xtype,'relational')
    DXX     = X;                        % Ignote the distance parameter            
else
    error('Unknown input data type. It should be "feature" or "relational".');
end
% ------------------------------------
% Evaluacion del IVG
DX = DXX.^2;
S = zeros(1,K);
for i = 1:K
    XC = DX(clust==cnames(i),clust==cnames(i));
    ai = sum(XC,2)/max(Nk(i)-1, 1); % modificacion 25-10-16
    %ai = mean(XC,2);
    dd = zeros(Nk(i),K-1);
    k = 0;
    for j = setdiff(1:K,i)
        k = k+1;
        XK = DX(clust==cnames(i),clust==cnames(j));
        dd(:,k) = mean(XK,2);
    end
    bi = min(dd,[],2);
    S(i) = sum((bi-ai)./max([ai bi],[],2),1);
end
f = (1/N)*sum(S);
