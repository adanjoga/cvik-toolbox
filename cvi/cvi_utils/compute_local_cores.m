
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2019 - Dongdong Cheng -- A Novel Cluster Validity Index Based on Local Cores
% LCCV index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [cores,short_path,local_core] = compute_local_cores(X,DXX)
% X might be used to print the local cores but it is not essential
[N,~]=size(X); 
[sdist,index]=sort(DXX,2); % Sort the distances

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NaN-Searching algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('Start running NaN-Searching algorithm...\n');
r=1;
flag=0;         
nb=zeros(1,N);  %The number of each point's reverse neighbor
count=0; count1=0;    
while flag==0
    for i=1:N
        k=index(i,r+1);
        nb(k)=nb(k)+1;
    end
    r=r+1;
    count2=0;
    for i=1:N
        if nb(i)==0
            count2=count2+1;
        end
    end
    if count1==count2
        count=count+1;
    else
        count=1;
    end
    if count2==0 || (r>2 && count>=2)   %The terminal condition
        flag=1;
    end
    count1=count2;
end

supk=r-1;               %The characteristic value of natural neighbor
max_nb=max(nb);         %The maximum value of nb
%fprintf('The characteristic value is %d\n',supk);
%fprintf('The maximum value of nb is %d\n',max_nb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate the density of each point
rho=zeros(N,1);
Non=max_nb;
for i=1:N
    d=0;
    for j=1:Non+1
        d=d+sdist(i,j);
    end
    rho(i)=(Non/d);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LORE Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,ordrho]=sort(rho,'descend');%sort the points according to the density
local_core = zeros(N,1);%The representative of each point
%fprintf('Starting running LORE algorithm...\n');
for i=1:N
         p=ordrho(i);
         maxrho=rho(p);
         maxindex=p;
         %Find the point with maximum density in the local neighbors
         for j=1:nb(p)+1
             x=index(p,j);
             if maxrho<rho(x)
                 maxrho=rho(x);
                 maxindex=x;
             end
         end
         %Assign representative of the point with maximum density
         if local_core(maxindex)==0
             local_core(maxindex)=maxindex;            
         end
         %Assign representative of the local neighbors
         for j=1:nb(p)+1
             if local_core(index(p,j))==0
                 local_core(index(p,j))=local_core(maxindex);
             else%Determine the representative according to RCR
                 q=local_core(index(p,j));
                 if DXX(index(p,j),q)>DXX(index(p,j),local_core(maxindex))% rho(q)<rho(local_core(maxindex))%
                     local_core(index(p,j))=local_core(maxindex);
                 end
             end 
             %Determine the representative according to RTR
             for m=1:N
                 if local_core(m)==index(p,j)
                     local_core(m)=local_core(index(p,j));
                 end
             end
         end
  
end
% Find the cores
cluster_number=0;
for i=1:N
    if local_core(i)==i
       cluster_number=cluster_number+1;
       cores(cluster_number)=i;
    end
end

%fprintf('The number of initial clusters is %d\n',cluster_number);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw the cores and initial cluster result
% plot(X(:,1),X(:,2),'.');
% hold on;
% for i=1:N
%     plot([X(i,1),X(local_core(i),1)],[X(i,2),X(local_core(i),2)]);
%     hold on;
% end
% %drawcluster2(X,cl,cluster_number+1);
% hold on;
% plot(X(local_core,1),X(local_core,2),'r*','MarkerSize',8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the graph
weight=zeros(N,N);
for i=1:N
    for j=2:supk+1
        x=index(i,j);
        weight(i,x)=DXX(i,x);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the shortest path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fprintf('Start computing the graph-based distance between local cores...\n');
short_path=zeros(cluster_number,cluster_number); %The shortest path between local cores
weight2=sparse(weight);
for i=1:cluster_number
     short_path(i,i)=0; 
    [D,~,~] = graphshortestpath(weight2,cores(i),'METHOD','Dijkstra'); 
%    [D,Z]=dijkstra2(weight,cores(i));
     for j=i+1:cluster_number
         short_path(i,j)=D(cores(j));
         if short_path(i,j)==inf
             short_path(i,j)=0;
         end
         short_path(j,i)=short_path(i,j);
     end
end
maxd=max(max(short_path));
for i=1:cluster_number
    for j=1:cluster_number
        if short_path(i,j)==0
            short_path(i,j)=maxd;
        end
    end
end
end