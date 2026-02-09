clc,clear
load ./glc_rmse.mat
% CSSPv3, COSIPY, CSSPv3_nofre

K = [[255 199 204];[145 217 253];[235 255 235]]/255;
figure(1)
set(gcf,'Units','centimeters','Position',[5,5,35.5,20.5],'color','white');
h = bar(rmse(:,1:2),'BarWidth',1);
set(h(1),'FaceColor',K(1,:));
set(h(2),'FaceColor',K(2,:));
set(h, 'LineWidth',1);
xlim([0.5 6.5])
ylim([0 5.5])
set(gca,'YTick',0:1:5,'fontsize',14,'FontName','Time New Roman')
glacier_name = {'Alaska', 'Central Europe', 'Asia South East', ...
    'Asia West', 'Himalayan Region', 'Iceland'};
set(gca, 'XTickLabel',glacier_name,'fontsize',18,'FontName','Times New Roman');
set(gca, 'XTickLabelRotation', 0);

xlabel('Glaciers','FontName','Time New Roman','fontsize',18)
ylabel('RMSE','FontName','Time New Roman','fontsize',18)
grid on
set(gca, 'GridLineStyle', '--', ...
    'GridColor', [0.7 0.7 0.7], ...
    'GridAlpha', 0.6, ...
    'MinorGridLineStyle', '--')

numgroups = size(rmse(:,1:2),1);
numbars = size(rmse(:,1:2),2);
groupwidth = min(0.8,numbars/(numbars+1.5));

offset_x = 0;
for i = 1:numbars
    for j = 1:numgroups
        y = rmse(j,i);
        x = j - groupwidth/2 + (2*i-1)*groupwidth/(2*numbars);
        if y > 0
            offset_y = 0.02;
            text(x + offset_x, y + offset_y, sprintf('%.2f', y), ...
                'FontSize', 16, 'FontName', 'Times New Roman', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom','Color','r','FontWeight','bold');
        else
            offset_y = -0.02;
            text(x + offset_x, y + offset_y, sprintf('%.2f', y), ...
                'FontSize', 16, 'FontName', 'Times New Roman', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top','Color','r','FontWeight','bold');
        end
    end
end

legend('CSSPv3','COSIPY','FontName','Time New Roman','Location','Northeast','fontsize',16)
n_title = title('Evaluation of Glacier Mass Loss','FontSize',20,'FontName','Time New Roman',...
    'HorizontalAlignment','left', 'Units', 'normalized');
set(n_title, 'Position', [0.005, 1.015, 0]);

set(gca, 'LineWidth', 1.5)
box off
ax = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
    'YAxisLocation','right','Color','none','XColor','k','YColor','k');  
set(ax,'XTick', [],'YTick', []); 
set(gca, 'LineWidth', 1.5)
print(figure(1),'./Fig13_GML','-dpng','-r300')


fprintf('%-20s | %-10s | %-10s | %-10s | %-15s | %-15s\n', ...
    'Glacier', 'CSSPv3', 'COSIPY', 'noFre', 'Δ vs COSIPY (%)', 'Δ vs noFre (%)');
imp_cosipy = (rmse(:,1) - rmse(:,2)) ./ rmse(:,2) * 100;
imp_nofre  = (rmse(:,1) - rmse(:,3)) ./ rmse(:,3) * 100;

for i = 1:size(rmse,1)
    fprintf('%-20s | %-10.2f | %-10.2f | %-10.2f | %-15.2f | %-15.2f\n', ...
        glacier_name{i}, ...   % glacier/station name
        rmse(i,1), ...         % CSSPv3 raw RMSE
        rmse(i,2), ...         % COSIPY raw RMSE
        rmse(i,3), ...         % noFre raw RMSE
        imp_cosipy(i), ...     % signed %
        imp_nofre(i));
end