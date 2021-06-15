
% Select a solution from pareto front
... (orthogonal distance method)

function [ind,coord,distance] = unsupervised(data)

[~,mx_x] = max(data(:,1));
[~,mx_y] = max(data(:,2));

po1 = data(mx_x,:);
po2 = data(mx_y,:);

x = [po1(1) po2(1)];
y = [po1(2) po2(2)];

m = (po2(2) - po1(2))/(po2(1) - po1(1));
b = y(1) - m*x(1);
m2 = -1/m;

x2 = zeros(numel(data(:,1)),2);
y2 = zeros(numel(data(:,1)),2);
distance = zeros(numel(data(:,1)),1);

for i = 1:numel(data(:,1))
    b2 = data(i,2) - m2*data(i,1); 
    A = [-m 1;-m2 1]; 
    vb = [b;b2];
    po3 = A\vb; 
    x2(i,:) = [data(i,1) po3(1)];
    y2(i,:) = [data(i,2) po3(2)];
    distance(i) = dist(data(i,:),po3); 
end

[~,ind] = max(distance);
coord.x = x2; 
coord.y = y2;
coord.po1 = po1; 
coord.po2 = po2; 

