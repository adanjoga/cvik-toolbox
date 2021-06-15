function plotclusters(X,T,Y,M)

d = size(X,2);

if 2 <= d && d <= 3
    figure('color',[1 1 1]);    
    kT = unique(T);    % agrupamiento verdadero
    kY = unique(Y);    % agrupamiento resultante
    tcolors = rand(max(numel(kT),numel(kY)),3);
    if d == 2
        hold on;
        for i = 1:numel(kY)
            idx = Y==kY(i);
            plot(X(idx,1),X(idx,2),'o','MarkerEdgeColor','k','MarkerFaceColor',tcolors(i,:),'MarkerSize',10);
        end
        plot(M(:,1),M(:,2),'h','MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'MarkerSize',20);
        hold off; box on;

    else
        hold on;
        for i = 1:numel(kY)
            idx = Y==kY(i);
            plot3(X(idx,1),X(idx,2),X(idx,3),'o','MarkerEdgeColor','k','MarkerFaceColor',tcolors(i,:),'MarkerSize',10);
        end
        plot3(M(:,1),M(:,2),M(:,3),'h','MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'MarkerSize',20);
        hold off; box on; view(3);
       
    end
else
   error('No se puede graficar'); 
end