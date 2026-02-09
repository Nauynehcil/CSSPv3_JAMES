clc,clear
load ('./CSSP_output.mat','*_cssp');
load ('./MsTMIP_output.mat','*_models');

kge = cat(3,kge_cssp,kge_models);
cc = cat(3,cc_cssp,cc_models);
rmse = cat(3,rmse_cssp,rmse_models);
bias = cat(3,bias_cssp,bias_models);

Heat_Color = [255,255,255; 235,235,235; 252,248,187; 254,191,136; 246,137,109; 226,085,105; 184,056,116; 131,037,127; 085,020,118; 044,016,076]/255;

KGE_UPPER = squeeze(prctile(reshape(kge,[],size(kge,3)),75,1))';
KGE_LOWER = squeeze(prctile(reshape(kge,[],size(kge,3)),25,1))';
KGE_MID = squeeze(prctile(reshape(kge,[],size(kge,3)),50,1))';

CC_UPPER = squeeze(prctile(reshape(cc,[],size(kge,3)),75,1))';
CC_LOWER = squeeze(prctile(reshape(cc,[],size(kge,3)),25,1))';
CC_MID = squeeze(prctile(reshape(cc,[],size(kge,3)),50,1))';

RMSE_UPPER = squeeze(prctile(reshape(rmse,[],size(kge,3)),75,1))';
RMSE_LOWER = squeeze(prctile(reshape(rmse,[],size(kge,3)),25,1))';
RMSE_MID = squeeze(prctile(reshape(rmse,[],size(kge,3)),50,1))';

BIAS_UPPER = squeeze(prctile(reshape(abs(bias),[],size(kge,3)),75,1))';
BIAS_LOWER = squeeze(prctile(reshape(abs(bias),[],size(kge,3)),25,1))';
BIAS_MID = squeeze(prctile(reshape(abs(bias),[],size(kge,3)),50,1))';

eva_mid = [KGE_MID,CC_MID,RMSE_MID,BIAS_MID];
eva_75 =  [KGE_UPPER,CC_UPPER,RMSE_UPPER,BIAS_UPPER];
eva_25 =  [KGE_LOWER,CC_LOWER,RMSE_LOWER,BIAS_LOWER];
i_eva_75 = eva_75 - eva_mid;
i_eva_25 = eva_mid - eva_25;

figure(1)
% 1 Africa 2 Asia 3 Europe 4 North America 5 South-West Pacific 6 South America
title_name = {'(a) KGE and CC of Different Models Compared with FLUXSat',...
              '(b) Absolute BIAS and RMSE of Different Models Compared with FLUXSat'};
K = [[255 199 204];[145 217 253];[235 255 235]]/255;
set(gcf,'Units','centimeters','Position',[5,5,44.5,28.5],'color','white')
na = tight_subplot(2,1,[.1 .045],[.055 .055],[.055 .055]);
axes(na(1))
h = bar(eva_mid(:,1:2),1); 
set(h(1),'FaceColor',K(1,:));
set(h(2),'FaceColor',K(2,:));
set(h, 'LineWidth',1);
set(gca,'YLim',[-0.3 1.2]);
yticks(-0.2:0.2:1)
xticks(1:10)
set(gca,'XTickLabel',{'CSSPv3','CLM4-VIC','CLM4','DLEM', ...
    'GTEC','ISAM','LPJ-WSL','ORCHIDEE-LSCE','SiBCASA','VEGAS2.1'},...
    'fontsize',15,'FontName','Times New Roman');
set(gca, 'XTickLabelRotation', 0); 
xlabel('Models','FontName','Time New Roman','fontsize',18);
ylabel('KGE & CC','FontName','Time New Roman','fontsize',18);
hold on
numgroups = size(eva_mid(:,1:2),1);
numbars =   size(eva_mid(:,1:2),2);
groupwidth = min(0.8,numbars/(numbars+1.5));
for i = 1:numbars
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth /(2*numbars);
    errorbar(x,eva_mid(:,i),i_eva_25(:,i),i_eva_75(:,i), 'k', 'linestyle', 'none', 'lineWidth', 1,'Color',K(i,:)-0.25);
end
offset_x = 0;
for i = 1:numbars
    for j = 1:numgroups
        y = eva_mid(j,i);
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
hold off
set(gca, 'LineWidth', 1)
legend('KGE','CC',...
       'FontName','Time New Roman','Location','Northeast','fontsize',14);
n_title = title(cell2mat(title_name(1)),'FontSize',18,'FontName','Time New Roman',...
    'HorizontalAlignment','left', 'Units', 'normalized');
set(n_title, 'Position', [0.005, 1.025, 0]); 
box off
ax = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
    'YAxisLocation','right','Color','none','XColor','k','YColor','k'); 
set(ax,'XTick', [],'YTick', []); 
set(gca, 'LineWidth', 1)

axes(na(2))
h = bar(eva_mid(:,3:4),1); 
set(h(1),'FaceColor',K(1,:));
set(h(2),'FaceColor',K(2,:));
set(h, 'LineWidth',1);
set(gca,'YLim',[0 80]);
yticks(0:20:80)
xticks(1:10)
set(gca,'XTickLabel',{'CSSPv3','CLM4-VIC','CLM4','DLEM', ...
    'GTEC','ISAM','LPJ-WSL','ORCHIDEE-LSCE','SiBCASA','VEGAS2.1'},...
    'fontsize',15,'FontName','Times New Roman');
set(gca, 'XTickLabelRotation', 0); 
xlabel('Models','FontName','Time New Roman','fontsize',18);
ylabel('RMSE & ABS(BIAS)','FontName','Time New Roman','fontsize',18);
hold on
numgroups = size(eva_mid(:,3:4),1);
numbars =   size(eva_mid(:,3:4),2);
groupwidth = min(0.8,numbars/(numbars+1.5));
for i = 1:numbars
    x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth /(2*numbars);
    errorbar(x,eva_mid(:,i+2),i_eva_25(:,i+2),i_eva_75(:,i+2), 'k', 'linestyle', 'none', 'lineWidth', 1,'Color',K(i,:)-0.25);
end
offset_x = 0;
for i = 1:numbars
    for j = 1:numgroups
        y = eva_mid(j,i+2);
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
hold off
set(gca, 'LineWidth', 1)
legend('ABS(BIAS)','RMSE',...
       'FontName','Time New Roman','Location','Northeast','fontsize',14);
n_title = title(cell2mat(title_name(2)),'FontSize',18,'FontName','Time New Roman',...
    'HorizontalAlignment','left', 'Units', 'normalized');
set(n_title, 'Position', [0.005, 1.025, 0]); 
box off
ax = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
    'YAxisLocation','right','Color','none','XColor','k','YColor','k');
set(ax,'XTick', [],'YTick', []);
set(gca, 'LineWidth', 1)

print(figure(1),'./Fig15_GPP_Com','-dpng','-r300')
