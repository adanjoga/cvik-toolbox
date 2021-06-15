function plotPFA(PFront, Nobj, idx1, idx2)

hold on; %axis square;
if Nobj == 2
    p1 = scatter(PFront(:,1), PFront(:,2),20,'o','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1.2); hold on;
    p1.MarkerFaceAlpha = 0.3; p1.MarkerEdgeAlpha = 0.4;
    xlabel('CVI 1'); ylabel('CVI 2');
    
    nPF = getParetoFrontLine(PFront);
    ppf = plot(nPF(:,1), nPF(:,2),'-k','LineWidth',1.2); hold on; ppf.Color(4) = 0.4;
    
    if nargin > 2 && (idx1 <= size(PFront,1))
        hold on;
        scatter(PFront(idx1,1), PFront(idx1,2),22,'o','MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1.2); hold on;
    end
    if nargin > 3 && (idx2 <= size(PFront,1))
        hold on;
        scatter(PFront(idx2,1), PFront(idx2,2),40,'o','MarkerEdgeColor','r','LineWidth',1.2); hold on;
    end    
    
elseif Nobj == 3
    
    p1 = scatter3(PFront(:,1), PFront(:,2), PFront(:,3), 30,'o','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',1.2); hold on;
    p1.MarkerFaceAlpha = 0.3; p1.MarkerEdgeAlpha = 0.4;
    
    if nargin > 2 && (idx1 <= size(PFront,1))
        hold on;
        scatter(PFront(idx1,1), PFront(idx1,2),PFront(idx1,3), 40,'o','MarkerFaceColor','b','MarkerEdgeColor','b','LineWidth',1.2); hold on;
    end
    if nargin > 3 && (idx2 <= size(PFront,1))
        hold on;
        scatter(PFront(idx2,1), PFront(idx2,2),PFront(idx1,3),60,'o','MarkerEdgeColor','r','LineWidth',1.2); hold on;
    end
end

if nargin == 3
    legend('PFA solution','Pareto front approx.','best solution')
end
if nargin == 4
    legend('PFA solution','Pareto front','best ARI', 'best Sil')
end

    
end

%% Function
function nPF = getParetoFrontLine(PF)    
    nPF = [];
    %firstp = [PF(1,1) PF(1,2)+1]; nPF = [nPF; firstp];
    for i=1:length(PF)-1
        nPF = [nPF; PF(i,:)];
        newp = [PF(i+1,1) PF(i,2)];
        nPF = [nPF; newp];
    end
    nPF = [nPF; PF(i+1,:)];
end   