clc,clear
load ./Metric_Data.mat
load ../mask.mat

M = [069,117,180; 224,243,248;
    255,245,229; 253,209,155; 252 140 089; 216,051,033; 135,000,000]/255;
N = [255,255,255; 235,235,235; 252,248,187; 254,191,136; 246,137,109; 226,085,105; 184,056,116; 131,037,127; 085,020,118; 044,016,076]/255;

CSSP_kge = kge(:,1);  % Extract KGE of CSSPv3 
figure(1)
set(gcf,'Units','centimeters','Position',[5,1.5,38.5,25],'color','white')
ha = tight_subplot(1,1,[.075 .075],[.505 .075],[.05 .05]);
ma = tight_subplot(1,2,[.075 .075],[.05 .55],[.105 .105]);
% Region names and corresponding lon/lat ranges
regions = {
    '(a) North America & West Europe', [-160, 40], [35, 75];
    '(b) East Asia',     [85, 145],   [20, 50];
    '(c) Australia',     [105, 165],  [-42, -8];
    };
for i = 1:3
    % Select subplot position
    if i == 1
       axes(ha(1))
    end
    if i >= 2
       axes(ma(i-1))
    end
    m_proj('Equidistant Cylindrical', 'long', regions{i,2}, 'lat', regions{i,3}); hold on
    m_coast('patch', [245, 250, 250]/255,'edgecolor', 'none');
    m_grid('xtick', floor(regions{i,2}(1)/10)*10 : 20 : ceil(regions{i,2}(2)/10)*10,...
           'ytick', floor(regions{i,3}(1)/10)*10 : 10 : ceil(regions{i,3}(2)/10)*10,...
           'tickdir','out','FontSize',12,'FontName','Times New Roman');
    for j = 1:size(Basin_SHP,1)
        m_line(Basin_SHP(j).Lon,Basin_SHP(j).Lat,'color',[0.5 0.5 0.5], 'linewidth', 1);
    end
    lon_range = regions{i,2};
    lat_range = regions{i,3};
    in_region = sta_lon >= lon_range(1) & sta_lon <= lon_range(2) & ...
                sta_lat >= lat_range(1) & sta_lat <= lat_range(2);
    m_scatter(sta_lon(in_region),sta_lat(in_region),30, CSSP_kge(in_region),'filled');
    colormap(M)
    clim([-0.4,1])

    n_title = title(regions{i,1},'FontSize',16,'FontName','Time New Roman',...
        'HorizontalAlignment','left', 'Units', 'normalized');
    set(n_title, 'Position', [0.005, 1.025, 0]);

    kge_median = median(CSSP_kge(in_region,1),'omitnan'); % Median Value of KGE
    text(1.0, 1.025,['Median: ',num2str(round(kge_median,2))], 'Units', 'normalized', ...
        'FontSize', 16, 'FontName', 'Times New Roman', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'Color', 'r','FontWeight','bold');

    if i == 3
        c = colorbar('Position', [0.92 0.1 0.015 0.8]);
        set(c,'ticks',(-0.2:0.2:0.8),'ticklength',0.015,'fontsize',14,'FontName','Times New Roman')
    end
end

model_names = {'CSSPv3','GLEAMv4.2','ERA5','ERA5\_LAND','GLDAS\_CLSM','GLDAS\_NOAH','GLDAS\_VIC',...
    'DBH','H08','LPJML','PCR-GLOBWB','WATERGAP','VIC','CLM4 (ISIMIP2a)',...
    'BIOME-BGC','CLASS-CTEM-N','CLM4-VIC','CLM4 (MsTMIP)','DLEM', ...
    'GTEC','ISAM','JPL-CENTURY','LPJ-WSL','ORCHIDEE-LSCE','SiBCASA','VEGAS2.1','VISIT'};
model_names = cellstr(model_names');

KGE_mid = round(squeeze(median(kge,1,'omitnan')),2);
RMSE_mid = round(squeeze(median(rmse,1,'omitnan')),2);
CC_mid = round(squeeze(median(cc,1,'omitnan')),2);

% Normalize metrics to [0,1] for heatmap comparison
KGE_norm = rescale(KGE_mid);
RMSE_norm = 1 - rescale(RMSE_mid);
CC_norm = rescale(CC_mid);
norm_data = [KGE_norm', CC_norm', RMSE_norm'];
norm_data_ori = [KGE_mid', CC_mid', RMSE_mid']; % for table sx

figure(2);
set(gcf, 'Units', 'centimeters', 'Position', [5,5,25,28], 'Color', 'white');
% Heatmap visualization of normalized performance
h = heatmap({'KGE', 'CC', 'RMSE'}, model_names, norm_data(:,1:3), ...
    'Colormap', N, ...
    'ColorLimits', [0 1], ....
    'CellLabelFormat', '%.2f');
h.FontName = 'Times New Roman';
h.FontSize = 13;
xlabel(h, 'Evaluation Metrics');
ylabel(h, 'Models');
% Add labels indicating best and worst performance
annotation(gcf,'textbox',[0.875,0.885,0.5,0.075],'String','Best','fontsize',12,'linestyle','none', 'Color', 'k');
annotation(gcf,'textbox',[0.875,0.025,0.5,0.075],'String','Worst','fontsize',12,'linestyle','none', 'Color', 'k');

print(figure(1),'./Fig9_FEVPA_KGE','-dpng','-r300')
print(figure(2),'./Fig11_FEVPA_Heatmap','-dpng','-r300')

