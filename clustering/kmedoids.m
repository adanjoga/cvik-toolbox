function cellout = kmedoids(xDist,k)
%kmedoids - an iterative algorithm to perform k-medoids partional clustering
%
% kmedoids implements an algorithm in which a k-means like update to cluster 
% assignments is done followed by an update of medoid given the current cluster
% assignments.

maxIterations = 200;
cellout = cell(8,1); % cellout{1} = total sum of distances
                     % cellout{2} = iteration
                     % cellout{3} = idx;
                     % cellout{4} = Medoids
                     
N = size(xDist,1);
medoidIndex = datasample(1:N,k,'Replace',false);

n = size(xDist,1);
for kter = 1:maxIterations
    nochange = false;
    for i = 1:k
        
        nomedoidIndex = setdiff(1:n,medoidIndex);
        [m,h] = min(xDist(:,medoidIndex),[],2);
        
        nd = xDist(:,nomedoidIndex);
        
        gain = sum(max(bsxfun(@minus,m(h~=i),nd(h~=i,:)),0),1);
        gaini = sum(m(h==i));
        
        %best cost of assigning outside samples to any other medoid
        mm = min(xDist(h==i,setdiff(medoidIndex,medoidIndex(i))),[],2);
        if ~isempty(mm)
            %gain -or lost- of reassigning inside samples to a new medoid
            gaini = gaini-sum(bsxfun(@min,mm,nd(h==i,:)),1);
        end
        %total gain
        gaint = gain+gaini;
        
        if any(gaint>0)
            [~,uu] = max(gaint);
            medoidIndex(i) = nomedoidIndex(uu);
            nochange = true;
        end
    end
    
    if ~nochange
        break
    end
    
end

DXM = xDist(:,medoidIndex);
[~,h] = min(DXM,[],2);

cellout{1} = sum(m);
cellout{2} = kter;
cellout{3} = h; % clusterID
cellout{4} = medoidIndex;

end

