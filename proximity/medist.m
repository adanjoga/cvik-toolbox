function GXX = medist(DXX)

% Construccion del grafo no dirigido (triangulo inferior de DXX)
LXX = tril(DXX);
UG = sparse(LXX);
% Encontrar el arbol de recubrimiento minimo
[MST,~] = graphminspantree(UG);

nVertex = length(MST);
MED = zeros(nVertex);
for i = 1:nVertex;
    for j = i+1:nVertex
        [~,path,~] = graphshortestpath(MST,i,j,'directed',false);
        v_i = path(1:end-1);    % Vertex i
        v_j = path(2:end);      % Vertex j
        subMST = MST(v_i,v_j) + MST(v_j,v_i);   % Subgraph of MST
        [~,~,distances] = find(subMST);   
        MED(i,j) = max(distances);  % MED distance from subgraph of MST
    end
end
GXX = MED + MED';