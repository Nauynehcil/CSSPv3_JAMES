clc,clear
load ./Metric_Data.mat
load ../mask.mat

% Given that many stations lack basin area information, particularly those in the China region,
% Basin categories are defined based on the long-term average monthly streamflow at each basin outlet.
% The KGE values are aggregated to assess model skill across different basin sizes.

% Perform streamflow classification (grouping by magnitude)
% Define thresholds based on the 33rd and 66th percentiles
edges = prctile(Mean_VAL, [33 66]); disp(edges);
% Categorize stations into groups (Small, Medium, Large flow)
groups = discretize(Mean_VAL, [-inf edges inf]);
groups_num = length(edges) + 1;
kge_groups = cell(1,3); 
area_group = cell(1,3);
mid_area_group = zeros(1,3);
for i = 1:groups_num
    % Extract KGE values and drainage areas for each flow group
    kge_groups{i} = kge(groups == i, :);
    area_group{i} = sta_area(groups == i, :);
    % Calculate the median area for each group to represent the scale
    mid_area_group (i) = median(area_group{i},'omitnan');
end
% Append the "All Stations" group as the final category
kge_groups{end+1} = kge; 
groups_num = groups_num + 1;

% Calculate KGE frequency distribution for each flow category
edges = [-Inf, 0, 0.25, 0.5, 0.75, Inf];
% Dimensions: [Groups (0-33%, 33-66%, 66-100%, All), Models (1-8), KGE Intervals]
perf_ratio = zeros(groups_num, 8, 5); 
median_kge = zeros(groups_num, 8); 
kgestr = cell(groups_num, 8); % Fixed index to match model count
for i = 1:groups_num
    % Calculate median KGE per group and model
    median_kge(i,:) = median(kge_groups{i}(:,:),1,'omitnan');
    for model = 1:8
        % Convert median KGE to string for labeling; handle negative values
        if median_kge(i,model) >= 0
            kgestr{i,model} = num2str(median_kge(i,model));
        else
            kgestr{i,model} = '<0';
        end
        % Calculate the percentage of stations within each KGE interval (Model Skill)
        vals = kge_groups{i}(:,model);
        counts = histcounts(vals, edges);
        perf_ratio(i, model, :) = counts / sum(counts);
    end
end
% 
figure(1)
basin_name = {'(a) All Watersheds',...
              '(b) Watersheds with Mean Streamflow between 0-33% Quantile',...
              '(c) Watersheds with Mean Streamflow between 33-66% Quantile',...
              '(d) Watersheds with Mean Streamflow between 66-100% Quantile'};
model_names = {'CSSPv3','CSSPv2','DBH','H08','LPJML','PCR-GLOBWB','WATERGAP','CLM4'};
set(gcf,'Units','centimeters','Position',[0.5,0.5,60.5,32.5],'color','white')
na = tight_subplot(2,2,[.085 .045],[.055 .115],[.045 .005]);

% Read ColorMap
Color_Var = [189,229,189; 254,218,185; 159,187,224; 254,132,196; 187,187,187]/255;
for i = 1:groups_num
    axes(na(i))
    if i == 1
        k = groups_num; % Plot All Station
    else
        k = i - 1; % Plot 0-33% -> 33-66% -> 66-100%
    end
    data = squeeze(perf_ratio(k,:,:)); % Groups
    b = bar(data,'stacked','BarWidth', 0.6);
    for j = 1:size(data,2)
        b(j).FaceColor = 'flat';
        b(j).CData = Color_Var(j,:);
        b(j).FaceAlpha = 0.7;
    end

    ylim([0 1.05])
    set(gca,'YTick',0:0.2:1,'fontsize',12,'FontName','Time New Roman')
    set(gca,'YTickLabel',{'0','20%','40%','60%','80%','100%'},'fontsize',12,'FontName','Time New Roman');

    xlim([0.5 7.5])
    xticks(1:8)
    set(gca,'XTickLabel',model_names,'fontsize',15,'FontName','Times New Roman');
    set(gca, 'XTickLabelRotation', 0);
    xlabel('Models','FontName','Time New Roman','fontsize',17);
    ylabel('Percentage(%)','FontName','Time New Roman','fontsize',17);
    n_title = title(cell2mat(basin_name(i)),'FontSize',18,'FontName','Time New Roman',...
        'HorizontalAlignment','left', 'Units', 'normalized');
    set(n_title, 'Position', [0.005, 1.025, 0]);
    box off
    ax = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
        'YAxisLocation','right','Color','none','XColor','k','YColor','k');
    set(ax,'XTick', [],'YTick', []);
    hold off
end

% Note: The 'Position' parameters below may require manual adjustment based on image dimensions and screen resolution.
% Create a horizontal legend for the KGE intervals
labels = {'<0','0-0.25','0.25-0.5','0.5-0.75','≥ 0.75'};
lgd = legend(b,labels, 'Orientation', 'horizontal', ...
    'FontSize', 17, 'FontName', 'Times New Roman', ...
    'Box', 'off', 'Units', 'normalized', ...
    'Position', [0.74,0.915,0.21,0.03]); 

% Add a text label 'KGE:' to the left of the legend
annotation('textbox', [0.665, 0.915, 0.06, 0.03], ...
    'String', 'KGE:', 'EdgeColor', 'none', ...
    'FontSize', 17, 'FontName', 'Times New Roman', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
print(figure(1),'./Fig4_Streamflow_Basins','-dpng','-r300')

% --- Print Median KGE by Category ---
group_labels = {'0-33%', '33-66%', '66-100%', 'Total'};
fprintf('\n[Statistical Summary: Median KGE]\n');
fprintf('%-10s |', 'Basin Class');
for m = 1:8
    fprintf(' %-10s |', model_names{m});
end
fprintf('\n%s\n', repmat('-', 1, 105)); % Compact separator
for i = 1:4
    fprintf('%-10s |', group_labels{i});
    for m = 1:8
        val = median_kge(i, m);
        if val < 0, fprintf(' %-10s |', '<0');
        else, fprintf(' %-10.2f |', val); end
    end
    fprintf('\n');
end

% --- Print KGE Interval Ratios (Percentage) ---
kge_intervals = {'<0', '0-0.25', '0.25-0.5', '0.5-0.75', '>=0.75'};
for i = 1:4
    fprintf('\n[KGE distribution over basins of: %s]\n', group_labels{i});
    fprintf('%-12s |', 'Model');
    for k = 1:5, fprintf(' %-10s |', kge_intervals{k}); end
    fprintf('\n%s\n', repmat('-', 1, 75));
    for m = 1:8
        fprintf('%-12s |', model_names{m});
        for k = 1:5
            val = perf_ratio(i, m, k) * 100;
            fprintf(' %-8.1f%% |', val);
        end
        fprintf('\n');
    end
end