clc,clear
load ./Representative_Sites.mat
load ../mask.mat

fig_value = cat(3,hum_obs,hum_sim);
figure(1)
letters = 'a':'z'; ti = arrayfun(@(x) ['(' x ')'], letters, 'UniformOutput', false);
set(gcf,'Units','centimeters','Position',[1,3,60,30.5],'color','white')
na = tight_subplot(4,2,[.055 .035],[.05 .05],[.045 .015]);
% Define Colormap
% Black: Observation | Red: CSSPv3 | Blue: CSSPv2
K = [000, 000, 000; 228, 026, 028;  055, 126, 184; 055, 228, 028;]/255 ;
x = 1:nmon;
for i = 1:8
    axes(na(i))
    plot_value = squeeze(fig_value(i,:,:))';
    % Identify valid data range (clip the plot to available observation period)
    valid = ~isnan(plot_value(1,:)); % if obs has value
    xs = find(valid,1,"first"); xe = find(valid,1,"last"); % Only plot where have value
    % Calculate dynamic Y-axis limits with a 25% overhead for text labels
    ys = min(plot_value,[],'all') * 0.95;  ye = max(plot_value,[],'all') * 1.25;
    for j = 3:-1:1
        p(j) = plot(x, plot_value(j,:),'Color', K(j,:), 'LineWidth',1.5); hold on
    end
    xlim([xs xe])
    xticks(12:24:240);
    ylim([ys ye])

    ax = gca; ax.YAxis.Exponent = 0;
    set(gca,'XTickLabel',1992:2:2010,'fontsize',12,'FontName','Times New Roman');
    grid on; set(gca,'GridLineStyle','--','GridColor',[0.7 0.7 0.7],'GridAlpha',0.6)

    ylabel('m^3/s/month')
    box off;
    ax = axes('Position',get(gca,'Position'),'XAxisLocation','top', ...
        'YAxisLocation','right','Color','none','XColor','k','YColor','k');
    set(ax,'XTick', [],'YTick', []);
    set(gca, 'FontSize', 13, 'LineWidth', 0.5, ...
        'XColor', 'k', 'YColor', 'k', 'ZColor', 'k')

    if i == 2
        legend(p, {'OBS','CSSPv3','CSSPv2'}, ... %
            'Location','northeast', ...
            'FontSize',14, ...
            'FontName','Times New Roman', ...
            'NumColumns',3);
    end
    n_title = title([cell2mat(ti(i)),' ',cell2mat(hum_sta_name(i))],'FontSize',16,'FontName','Time New Roman',...
        'HorizontalAlignment','left', 'Units', 'normalized');
    set(n_title, 'Position', [0.005, 1.025, 0]);

    text(1.0, 1.025,hum_basin_name(i), 'Units', 'normalized', ...
        'FontSize', 16, 'FontName', 'Times New Roman', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

    text(0.01, 0.9,['CSSPv3: KGE = ',num2str(hum_kge(i,1),2)], 'Units', 'normalized', ...
        'FontSize', 14, 'FontName', 'Times New Roman', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','Color',K(2,:));

    text(0.01, 0.8,['CSSPv2: KGE = ',num2str(hum_kge(i,2))], 'Units', 'normalized', ...
        'FontSize', 14, 'FontName', 'Times New Roman', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','Color',K(3,:));
end
print(figure(1),'./Fig5_Streamflow_Hum','-dpng','-r300')
clearvars fig_value

fig_value = cat(3,glc_obs,glc_sim);
figure(2)
letters = 'a':'z'; ti = arrayfun(@(x) ['(' x ')'], letters, 'UniformOutput', false);
set(gcf,'Units','centimeters','Position',[1,3,60,24.5],'color','white')
na = tight_subplot(3,2,[.075 .035],[.05 .05],[.045 .015]);
% Define Colormap
% Black: Observation | Red: CSSPv3 | Blue: CSSPv2
K = [000, 000, 000; 228, 026, 028;  055, 126, 184; 055, 228, 028;]/255 ;
x = 1:nmon;
for i = 1:6
    axes(na(i))
    plot_value = squeeze(fig_value(i,:,:))';
    % Identify valid data range (clip the plot to available observation period)
    valid = ~isnan(plot_value(1,:)); % if obs has value
    xs = find(valid,1,"first"); xe = find(valid,1,"last"); % Only plot where have value
    % Calculate dynamic Y-axis limits with a 25% overhead for text labels
    ys = min(plot_value,[],'all') * 0.95;  ye = max(plot_value,[],'all') * 1.25;
    for j = 3:-1:1
        p(j) = plot(x, plot_value(j,:),'Color', K(j,:), 'LineWidth',1.5); hold on
    end
    xlim([xs xe])
    xticks(12:24:240);
    ylim([ys ye])

    ax = gca; ax.YAxis.Exponent = 0;
    set(gca,'XTickLabel',1992:2:2010,'fontsize',12,'FontName','Times New Roman');
    grid on; set(gca,'GridLineStyle','--','GridColor',[0.7 0.7 0.7],'GridAlpha',0.6)

    ylabel('m^3/s/month')
    box off;
    ax = axes('Position',get(gca,'Position'),'XAxisLocation','top', ...
        'YAxisLocation','right','Color','none','XColor','k','YColor','k');
    set(ax,'XTick', [],'YTick', []);
    set(gca, 'FontSize', 13, 'LineWidth', 0.5, ...
        'XColor', 'k', 'YColor', 'k', 'ZColor', 'k')

    if i == 2
        legend(p, {'OBS','CSSPv3','CSSPv2'}, ... %
            'Location','northeast', ...
            'FontSize',14, ...
            'FontName','Times New Roman', ...
            'NumColumns',3);
    end
    n_title = title([cell2mat(ti(i)),' ',cell2mat(glc_sta_name(i))],'FontSize',16,'FontName','Time New Roman',...
        'HorizontalAlignment','left', 'Units', 'normalized');
    set(n_title, 'Position', [0.005, 1.025, 0]);

    text(1.0, 1.025,glc_basin_name(i), 'Units', 'normalized', ...
        'FontSize', 16, 'FontName', 'Times New Roman', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

    text(0.01, 0.9,['CSSPv3: KGE = ',num2str(glc_kge(i,1),2)], 'Units', 'normalized', ...
        'FontSize', 14, 'FontName', 'Times New Roman', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','Color',K(2,:));

    text(0.01, 0.8,['CSSPv2: KGE = ',num2str(glc_kge(i,2))], 'Units', 'normalized', ...
        'FontSize', 14, 'FontName', 'Times New Roman', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','Color',K(3,:));
end
print(figure(2),'./Fig6_Streamflow_Glc','-dpng','-r300')
