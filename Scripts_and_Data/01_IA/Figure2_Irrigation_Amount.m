clc,clear
load ../mask.mat
load ./IRRGDMD.mat

model_names = {'cssp','h08','lpjml','pcrglobwb','watergap','clm5'};
Pbias_edges = [-inf,-75,-25,25,75,inf];
counts_pct = zeros(length(Pbias_edges)-1, numel(model_names));
pbias = struct();
eva = zeros(2,numel(model_names));

for i = 1:numel(model_names)
    m = model_names{i};
    eval([m '(isnan(obs)) = nan;']);
    % Calculate Pbias
    pbias.(m) = 100 * (eval(m) - obs) ./ obs;
    % Calculate frequency distribution of Pbias
    flag = discretize(pbias.(m), Pbias_edges, 'IncludedEdge', 'right');
    valid = ~isnan(flag);
    counts_pct(:,i) = accumarray(flag(valid), 1, ...
        [length(Pbias_edges)-1, 1]) / length(clm5);
    % Evaluation metrics
    eva(1,i) = RMSE(eval(m), obs);
    eva(2,i) = BIAS(eval(m), obs);
end

% Performance improvement of CSSPv3 relative to other models
eva_up = (abs(eva(:,1)) - abs(eva)) ./ abs(eva) * 100;
% Remove grids with zero irrigation amount from plotting
xcssp(xcssp == 0) = nan;
% Redefined class edges for irrigation amount visualization
edges_map = [0.01,0.02,0.05,0.08,0.2,0.4,0.6,0.8];
% Reclassification for plotting (function defined at the end)
fig_map = reclass(xcssp,edges_map);

figure(1)
set(gcf,'Units','centimeters','Position',[5,5,30.5,27.5],'color','white');
ha = tight_subplot(1,1,[.1 .1],[.525 .075],[.075 .025]);
ma = tight_subplot(1,2,[.1 .1],[.075 .555],[.075 .075]);
box off

title_name = {'(a) CSSPv3-Simulated Irrigation Amount', ...
              '(b) The Comparison Between CSSPv3 and OBS', ...
              '(c) Percentage of Bias (PBIAS)'};
% Define colormaps
M = [255,245,229; 252,248,187; 254,191,136; 246,137,109; ...
    226,085,105; 184,056,116; 131,037,127; 085,020,118; 044,016,076]/255;
K = [254,132,196; 159,187,224; 189,229,189; ...
    254,218,185; 222,202,229; 187,187,187]/255;

% Subplot (a): spatial distribution of irrigation amount
axes(ha(1))
m_proj('Equidistant Cylindrical','long',[-180 180],'lat',[-60 90]); hold on
m_pcolor(xlon,xlat,fig_map); hold on
m_grid('xtick',(-180:30:180),'ytick',(-60:30:90), ...
    'tickdir','out','FontSize',12,'FontName','Times New Roman');
colormap(M)
% Plot irrigation district boundaries
for j = 1:size(Irr_SHP,1)
    m_line(Irr_SHP(j).Lon,Irr_SHP(j).Lat,'color',[0.5 0.5 0.5], 'linewidth', 1);
end
h_title = title(cell2mat(title_name(1)),'FontSize',15,'FontName','Time New Roman',...
    'HorizontalAlignment','left','Units', 'normalized');
set(h_title, 'Position',[0.005, 1.025, 0]);
clim([0,length(edges_map)+1])
c = colorbar;
set(c,'ticks',(0:1:length(edges_map)+1), ...
    'ticklength',0.015,'fontsize',13,'FontName','Times New Roman')
set(c,'ticklabels',{'0',edges_map,''}, ...
    'ticklength',0.015,'fontsize',13,'FontName','Times New Roman')
text(1, 1.025,'Unit: km^3/yr', 'Units', 'normalized', ...
    'FontSize', 15, 'FontName', 'Times New Roman', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

% Subplot (b): CSSPv3 vs OBS scatter plot
axes(ma(1))
xx = obs; yy = cssp;
valid = ~isnan(xx) & ~isnan(yy);
p = polyfit(xx(valid), yy(valid), 1);
slope = p(1);
cc = CC(yy,xx);
rmse = RMSE(yy,xx);
bias = BIAS(yy,xx);
bias_pct = PBIAS(yy,xx);
X = [ones(size(yy)) xx];
[b,bint,r,rint,stats] = regress(yy,X);
m = [(0:0.0001:1),1:1:10000];
n = m;
scatter(xx, yy, 20, 'o','filled'); hold on
plot(m,n,'k--', 'LineWidth', 1.5);
set(gca,'XTick',[0 0.0001 0.001 0.01 0.1 1 10 100 1000 10000], ...
    'YTick',[0 0.0001 0.001 0.01 0.1 1 10 100 1000 10000])
set(gca,'XLim',[0.0001 10000],'YLim',[0.0001 10000])
set(gca,'XScale','log','YScale','log','FontSize',13)
set(gca,'Fontname','Times New Roman','Box','on')
xlabel('Observed Irrigation Amount (km^3/yr)','FontSize',14)
ylabel('Simulated Irrigation Amount (km^3/yr)','FontSize',14)
m_title = title(cell2mat(title_name(2)), ...
    'FontSize',15,'FontName','Time New Roman', ...
    'HorizontalAlignment','left','Units', 'normalized');
set(m_title, 'Position',[0.005, 1.025, 0]);
text(0.4, 0.16,['RMSE = ',num2str(round(rmse,2)),' km^3/yr'], 'Units', 'normalized', ...
    'FontSize', 13, 'FontName', 'Times New Roman', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
text(0.4, 0.10,['Bias = ',num2str(round(bias,2)),' km^3/yr (',num2str(round(bias_pct,2)),'%)'], 'Units', 'normalized', ...
    'FontSize', 13, 'FontName', 'Times New Roman', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
text(0.4, 0.04,['Fitting: Slope = ',num2str(round(slope,3)),', R^2 = ',num2str(round(stats(1),2))], 'Units', 'normalized', ...
    'FontSize', 13, 'FontName', 'Times New Roman', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
text(0.8, 0.875,'1:1', 'Units', 'normalized', ...
    'FontSize', 13, 'FontName', 'Times New Roman', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

% Subplot (c): PBIAS distribution
axes(ma(2))
models = {'CSSPv3','H08','LPJML','PCR-GLOBWB','WATERGAP','CLM5'};
interval_labels = {'< -75%','-75% ~ -25%','-25% ~ 25%', ...
    '25% ~ 75%','> 75%'};
h = bar(counts_pct,'BarWidth',1); hold on
% Set bar colors
for i = 1:numel(h)
    h(i).FaceColor = K(i,:);
    h(i).EdgeColor = 'k';
    h(i).LineWidth = 1;
end
set(gca,'XTickLabel',interval_labels)
xlabel('Pbias Intervals')
ylabel('Probability Distribution')
xlim([0.5 5.5])
ylim([0 max(counts_pct(:))*1.2])
set(gca,'YTick',0:0.1:1,'fontsize',12,'FontName','Time New Roman')
set(gca,'YTickLabel',{'0','10%','20%','30%','40%','50%','60%','80%','100%'},'fontsize',11,'FontName','Time New Roman');
grid on
set(gca, 'GridLineStyle', '--', ...
    'GridColor', [0.7 0.7 0.7], ...
    'GridAlpha', 0.6, ...
    'MinorGridLineStyle', '--');
legend(models)
box on
set(gca, 'FontSize', 13, 'LineWidth', 1.2, ...
    'XColor', 'k', 'YColor', 'k', 'ZColor', 'k')
n_title = title(cell2mat(title_name(3)),'FontSize',15,'FontName','Time New Roman',...
    'HorizontalAlignment','left','Units', 'normalized');
set(n_title, 'Position',[0.005, 1.025, 0]);
print(figure(1),'./Fig2_IrrDmd','-dpng','-r300')

fprintf('%-15s | %-10s | %-10s | %-15s | %-15s\n', ...
       'Model Name', 'RMSE', 'BIAS', 'RMSE Improve (%)', 'BIAS Improve (%)');
model_names = {'CSSPv3','H08','LPJML','PCR-GLOBWB','WATERGAP','CLM5'};
for i = 1:6
    if i == 1
        rmse_imp_str = '-';
         bias_imp_str = '-';
    else
        rmse_imp = (eva(1,i) - eva(1,1)) / eva(1,i) * 100;
        rmse_imp_str = [num2str(round(rmse_imp, 1)), '%'];
        bias_imp = (abs(eva(2,i)) - abs(eva(2,1))) / abs(eva(2,i)) * 100;
        bias_imp_str = [num2str(round(bias_imp, 1)), '%'];
    end
    fprintf('%-15s | %-10.2f | %-10.2f | %-15s | %-15s\n', ...
        model_names{i}, eva(1,i), eva(2,i), rmse_imp_str, bias_imp_str);
end

function output = reclass(input,edges)
% Reclassify continuous values into discrete classes for plotting

output = input;
output(output == 0) = nan;

% First class: input <= edges(1)
mask = (input <= edges(1));
output(mask) = 0.5;

% Middle classes: edges(i-1) < input <= edges(i)
for i = 2:length(edges)
    val = i - 0.5;
    low = edges(i-1);
    high = edges(i);
    mask = (input > low & input <= high);
    output(mask) = val;
end

% Last class: input > edges(end)
output(input > edges(end)) = length(edges) + 0.5;
end
