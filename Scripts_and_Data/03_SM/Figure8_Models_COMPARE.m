clc,clear
load ./Metric_Data.mat
load ../mask.mat

for i = 1:7 % Different Region
    if i < 7
        pos = find(sta_region == i);
    else
        pos = find(sta_region > 0);
    end
    for j = 1:3 % Different Layer
        % Count number of stations with valid metrics for the current region/layer
        counts(i,j) = max(sum(~isnan(kge(pos,j,:)),1));
        % Calculate Median, 75th percentile, and 25th percentile (Interquartile Range)
        eva_mid(i,j,:) = median(kge(pos,j,:),'omitnan');
        eva_75(i,j,:) = prctile(kge(pos,j,:),75,1);
        eva_25(i,j,:) = prctile(kge(pos,j,:),25,1);
    end
end
% Calculate upper and lower error bar lengths relative to the median
i_eva_75 = eva_75-eva_mid;
i_eva_25 = eva_mid-eva_25;

% Sort regions based on the total number of stations for plotting priority
cnum = sum(counts,2);
[~,sort_idx] = sort(cnum,'descend');
region_names = {'Africa', 'Asia', 'Europe', 'North America', 'South-West Pacific', 'South America', 'Global'};
sorted_regions = region_names(sort_idx);

figure(1)
basin_name = {'(a) Global Stations','(b) Asia Stations',...
              '(c) North America Stations','(d) Europe Stations',...
              '(e) South-West Pacific Stations','(f) Africa Stations'};
K = [[255 199 204];[145 217 253];[235 255 235]]/255;
set(gcf,'Units','centimeters','Position',[5,5,62.5,34.5],'color','white')
na = tight_subplot(3,2,[.085 .045],[.055 .055],[.045 .005]);
for region = 1:6
    axes(na(region))

    % Extract median and error values for the current region
    plot_mid = squeeze(eva_mid(sort_idx(region),1:3,:))';
    plot_upper = squeeze(i_eva_75(sort_idx(region),1:3,:))';
    plot_bottom = squeeze(i_eva_25(sort_idx(region),1:3,:))';

    h = bar(plot_mid,1);
    set(h(1),'FaceColor',K(1,:));
    set(h(2),'FaceColor',K(2,:));
    set(h(3),'FaceColor',K(3,:));
    set(h, 'LineWidth',1);
    set(gca,'YLim',[-0.5 1]);
    set(gca,'YTick',-0.5:0.25:1,'fontsize',14,'FontName','Time New Roman')
    xticks(1:7)
    set(gca,'XTickLabel',{'CSSPv3','ESA-CCIv07.1','ERA5','ERA5\_LAND','GLDAS\_CLSM','GLDAS\_NOAH',...
        'GLDAS\_VIC'},'fontsize',14,'FontName','Times New Roman');
    set(gca, 'XTickLabelRotation', 0);
    xlabel('Models','FontName','Time New Roman','fontsize',16);
    ylabel('KGE','FontName','Time New Roman','fontsize',16);
    hold on

    numgroups = size(plot_mid,1);
    numbars = size(plot_mid,2);
    groupwidth = min(0.8,numbars/(numbars+1.5));
    for i = 1:numbars
        x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth /(2*numbars);
        errorbar(x,plot_mid(:,i),plot_bottom(:,i),plot_upper(:,i), 'k', 'linestyle', 'none', 'lineWidth', 1,'Color',K(i,:)-0.25);
    end
    offset_x = 0;

    % Intelligent Data Labeling
    % Logic to prevent text overlap when median values are too close
    for j = 1:numgroups
        y = plot_mid(j,:);
        offset_y = 0.125 * (y >= 0) + (-0.05) * (y < 0);
        thres = 0.075; % Proximity threshold for overlap detection
        for k = 1:2
            if abs(y(k) - y(k+1)) < thres
                if y(k) >= 0 && y(k+1) >= 0
                    % Both values are positive: bump up the larger one
                    if y(k) > y(k+1)
                        offset_y(k) = offset_y(k) + 0.045;
                    else
                        offset_y(k+1) = offset_y(k+1) + 0.045;
                    end
                elseif y(k) < 0 && y(k+1) < 0
                    % Both values are negative: bump down the smaller one
                    if y(k) < y(k+1)
                        offset_y(k) = offset_y(k) - 0.025;
                    else
                        offset_y(k+1) = offset_y(k+1) - 0.025;
                    end
                end
            end
        end
        % Place the formatted KGE text labels
        for i = 1:numbars
            x = j - groupwidth/2 + (2*i-1)*groupwidth/(2*numbars);
            text(x + offset_x, y(i) + offset_y(i) , sprintf('%.2f', y(i)), ...
                'FontSize', 15, 'FontName', 'Times New Roman', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top','Color','r','FontWeight','bold');
        end
    end
    hold off
    set(gca, 'LineWidth', 1)
    if region == 2
        legend('Surface Layer','Middle Layer','Deep Layer',...
            'FontName','Time New Roman','Location','Northeast','fontsize',16);
    end
    n_title = title(cell2mat(basin_name(region)),'FontSize',20,'FontName','Time New Roman',...
        'HorizontalAlignment','left', 'Units', 'normalized');
    set(n_title, 'Position', [0.005, 1.025, 0]);
    grid on
    set(gca, 'GridLineStyle', '--', ...
        'GridColor', [0.7 0.7 0.7], ...
        'GridAlpha', 0.6, ...
        'MinorGridLineStyle', '--')
    box off
    ax = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
        'YAxisLocation','right','Color','none','XColor','k','YColor','k'); 
    set(ax,'XTick', [],'YTick', []); 
    set(gca, 'LineWidth', 1)
end

print(figure(1),'./Fig8_SM_Models_KGE_Region','-dpng','-r300')


% Calculate summary statistics (Median) for KGE and RMSE across all stations
mid_kge = squeeze(median(kge,1,'omitnan'));
mid_rmse = squeeze(median(rmse,1,'omitnan'));
mid_rv = squeeze(median(rv,1,'omitnan'));

depth_labels = {'Surface', 'Middle', 'Root'};
model_names = {'CSSPv3','ESA-CCIv07.1','ERA5','ERA5\_LAND','GLDAS\_CLSM','GLDAS\_NOAH','GLDAS\_VIC'};
for j = 1:3
    fprintf('%-15s | %-15s | %-15s | %-15s | %-15s | %-18s | %-15s \n', ...
            cell2mat(depth_labels(j)), 'Median KGE', 'Median RMSE', 'Median RV', 'KGE Increase', 'RMSE Improve (%)', 'Relative Variability');
    for i = 1:7
        if i == 1
            kge_inc_str = '-';
            rmse_imp_str = '-';
        else
            kge_inc(i) = mid_kge(j,1) - mid_kge(j,i);
            kge_inc_str = num2str(round(kge_inc(i), 2));
            
            rmse_imp(i) = (mid_rmse(j,i) - mid_rmse(j,1)) / mid_rmse(j,i) * 100;
            rmse_imp_str = [num2str(round(rmse_imp(i), 1)), '%'];
        end
        fprintf('%-15s | %-15.2f | %-15.3f | %-15.2f | %-15s | %-18s |%-15.2f \n' , ...
            model_names{i}, mid_kge(j,i), mid_rmse(j,i), mid_rv(j,i), kge_inc_str, rmse_imp_str, mid_rv(j,i));
    end
    fprintf('\n') 
end

