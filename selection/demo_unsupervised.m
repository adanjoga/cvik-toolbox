clearvars; close all; clc;

load pareto_front.mat;

[ind,coord,distance] = unsupervised(data);

f1 = data(ind,1);
f2 = data(ind,2);

f = figure; 
scatter(data(:,1),data(:,2),50,'o',...
     'MarkerFaceColor',[1 0.1 0.3],...
     'MarkerEdgeColor',[1 0.1 0.3],...
     'MarkerEdgeAlpha',0.5,...
     'MarkerFaceAlpha',0.5);
xlabel('\textit{Hellinger distance} ($HD$)',...
       'Interpreter','LaTex',...
       'FontSize',13); 
ylabel('\textit{Cut-off point ($\bar{z}$)}',...
       'Interpreter','LaTex',...
       'FontSize',13);

axis([0 1 0 1]);
axis square;
po1 = coord.po1; 
po2 = coord.po2; 
x2 = coord.x;
y2 = coord.y; 
x = [po1(1) po2(1)];
y = [po1(2) po2(2)];
line(x,y,'Color','black','LineStyle','--','Linewidth',1.25);
for i = 1:numel(data(:,1))
    hold on;
    line(x2(i,:),y2(i,:),'Color',[0.5,0.5,0.5],'LineStyle',':','Linewidth',0.25); 
end
line(x2(ind,:),y2(ind,:),'Color',[0.5 0.5 0.5],'LineStyle','--','Linewidth',2.5);
hold on; 
scatter(x2(ind,1),y2(ind,1),400,'p',...
        'MarkerFaceColor',[1 0.1 0.3],...
        'MarkerEdgeColor',[1 0.1 0.3],...
        'MarkerEdgeAlpha',1,...
        'MarkerFaceAlpha',1);
box on; 
grid minor;
f.Position = [800 200 500 500];
f.Children.Position = [0.075 0.075 0.875 0.875];
f.Children.TickLength = [0.0075 0.0075];
[xlim1] = f.Children.XLim;
[ylim1] = f.Children.YLim;
f.Children.XTick = linspace(0,xlim1(2),10);
f.Children.YTick = linspace(0,ylim1(2),10);

legend1 = strcat('\textit{Pareto front}');
legend2 = strcat('\textit{Selected solution}');
c = get(f.Children,'Children');
[leg,icons] = legend(c([end 1]),...
                     legend1,legend2,...
                     'interpreter','latex');
leg.Title.Interpreter = 'latex';
leg.FontSize = 14;
icons(4).Children.FaceAlpha = 0.1; 
icons(4).Children.MarkerSize = 11;
drawnow; pause(1)
xlim([0 1]);