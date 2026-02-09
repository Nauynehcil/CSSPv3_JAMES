clc,clear
load ./Representative_Sites.mat
load ../mask.mat

fig_value = cat(3,target_mon,target_sim);
figure(1)
letters = 'a':'z'; ti = arrayfun(@(x) ['(' x ')'], letters, 'UniformOutput', false);
set(gcf,'Units','centimeters','Position',[1,3,40,25.5],'color','white')
na = tight_subplot(3,1,[.065 .035],[.05 .05],[.045 .015]); 
% Define Colormap
% Black: Observation | Red: CSSPv3 | Blue: CSSPv2
K = [000, 000, 000; 228, 026, 028; 055, 126, 184;]/255 ;
x = 1:288;
for i = 1:3
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
    xticks(1:12:288);
    ylim([ys ye])

    ax = gca; ax.YAxis.Exponent = 0;
    set(gca,'XTickLabel',1991:1:2020,'fontsize',12,'FontName','Times New Roman');
    grid on; set(gca,'GridLineStyle','--','GridColor',[0.7 0.7 0.7],'GridAlpha',0.6)

    ylabel('mm/month')
    box off;
    ax = axes('Position',get(gca,'Position'),'XAxisLocation','top', ...
        'YAxisLocation','right','Color','none','XColor','k','YColor','k');
    set(ax,'XTick', [],'YTick', []);
    set(gca, 'FontSize', 13, 'LineWidth', 0.5, ...
        'XColor', 'k', 'YColor', 'k', 'ZColor', 'k')

    if i == 1
        legend(p, {'OBS','CSSPv3','CSSPv2'}, ... %
            'Location','northeast', ...
            'FontSize',14, ...
            'FontName','Times New Roman', ...
            'NumColumns',3);
    end
    n_title = title([cell2mat(ti(i)),' ',cell2mat(target_name(i))],'FontSize',16,'FontName','Time New Roman',...
        'HorizontalAlignment','left', 'Units', 'normalized');
    set(n_title, 'Position', [0.005, 1.025, 0]);

    text(0.01, 0.9,['CSSPv3: KGE = ',num2str(target_kge(i,1),2)], 'Units', 'normalized', ...
        'FontSize', 14, 'FontName', 'Times New Roman', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','Color',K(2,:));

    text(0.01, 0.8,['CSSPv2: KGE = ',num2str(target_kge(i,2),2)], 'Units', 'normalized', ...
        'FontSize', 14, 'FontName', 'Times New Roman', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom','Color',K(3,:));

    text(1.0, 1.025,'ET', 'Units', 'normalized', ...
        'FontSize', 16, 'FontName', 'Times New Roman', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
end

print(figure(1),'./Fig10_FEVPA_IRRG','-dpng','-r300')