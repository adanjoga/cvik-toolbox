% -------------------------------------------------------------------------
%Aim: The matlab code of "An internal validity index based on density-involved distance"
%compute the Density-involved distance DD
% -------------------------------------------------------------------------
%Input:
%DXX: the Euclidean distance between objects in X
%Knn: the number of neighborhoods 
% -------------------------------------------------------------------------
%Output:
%results: the DD of d
% -------------------------------------------------------------------------
% Written by Lianyu Hu
% Department of Computer Science, Ningbo University 
% August 2018

function DD = density_involved_distance(DXX, Knn)
    %% compute density-involved distance
    N = length(DXX);
    % Knn = 7;
    [KNNG]=KNearestNeighborGraph(DXX,Knn);
    Den = zeros(N,1);
    for j = 1:N
        Den(j,1) = sum(DXX(j,KNNG{j,1}))/Knn;
    end
    fDen = Den/max(Den); %absolute-density distance factor
    Rel = repmat(Den,1,N)./repmat(Den',N,1);
    tmp1 = Rel + Rel';
    fRel = 1-exp(-abs(tmp1-2)); %relative-density distance factor
    nD = repmat(Den,1,N) + repmat(Den',N,1);
    relD = nD.*fRel;
    drD = DXX + relD; %directly density-reachable distance
    conD =  fast_PathbasedDist(drD); %connectivity distance
    tmp2 = sqrt(repmat(fDen,1,N).*repmat(fDen',N,1));
    DD = conD.*tmp2; %density-involved distance
end

