% u = clustering labeling
% v = true labeling
function out = inftheoryindex(u,v)
u = u(:); v = v(:);
n = numel(u);
uk = max(u);
vk = max(v);
cm = full(sparse(u,v,1,uk,vk));
rt = sum(cm,2);
ct = sum(cm,1);
H1 = -sum((rt/n).*log2((rt/n))); % Entropia de A
H2 = -sum((ct/n).*log2((ct/n))); % Entropia de B
% rtct = (rt/n)*(ct/n);
nz   = (cm/n)+(cm==0);
H12 = -sum(sum((cm/n).*log2(nz))); % Entropia conjunta
MI = H1+H2-H12; %sum(sum((cm./n).*log2(nz./rtct))) % Informacion mutua
NMI = (2*MI)/(H1+H2);   % Informacion mutua normalizada
VI = H1+H2-2*MI;    % Variation of information
NVI = VI/log2(n);   % Normalized variation of information
% Salidas
%out = [MI,NMI,VI,NVI];
out = NMI;