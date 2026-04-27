clc,clear
load Metric_Data.mat
load ../mask.mat

% Custom colormap for KGE: transitioning from Blue (Low/Poor) to Dark Red (High/Good)
Color_KGE = [069,117,180; 224,243,248; 255,245,229; 253,209,155; 252 140 089; 216,051,033; 135,000,000]/255;

% title name
title_name = {
    '(a) KGE of CSSPv3 Simulated Streamflow',...
    '(b) KGE of CSSPv2 Simulated Streamflow',...
    '(c) KGE of DBH Simulated Streamflow',...
    '(d) KGE of H08 Simulated Streamflow',...
    '(e) KGE of LPJML Simulated Streamflow',...
    '(f) KGE of PCR-GLOBWB Simulated Streamflow',...
    '(g) KGE of WATERGAP Simulated Streamflow',...
    '(h) KGE of CLM4 Simulated Streamflow'};

% Binning Parameters
res = 1; % Grid resolution (degrees)
lon_edges = -180:res:180;
lat_edges = -56:res:84;

figure(1)
set(gcf,'Units','centimeters','Position',[1.5,1.5,38.5,36],'color','white')

ha = tight_subplot(4,2,[.05 .035],[.055 .055], [.035 .005]);
for i = 1:8
    axes(ha(i))
    % Map Projection
    m_proj('Equidistant Cylindrical','long',[-180 180],'lat',[-60 90]); hold on
    m_grid('xtick',(-180:30:180),'ytick',(-60:30:90),'tickdir','out','FontSize',12,'FontName','Times New Roman');
    m_coast('patch',[250,250,250]/255,'edgecolor', 'none');

    % Spatial Binning Logic
    % Map station coordinates to grid indices
    binX = discretize(sta_lon, lon_edges);
    binY = discretize(sta_lat, lat_edges);
    sz = [length(lat_edges)-1, length(lon_edges)-1];
    valid = ~isnan(binX) & ~isnan(binY);

    % Calculate median KGE per Binning cell
    binned_kge = accumarray([binY(valid), binX(valid)], kge(valid, i), sz, @(x) median(x, 'omitnan'), NaN);
    [Xg, Yg] = meshgrid(lon_edges(1:end-1) + res/2, lat_edges(1:end-1) + res/2);

    % Plotting
    m_pcolor(Xg, Yg, binned_kge);
    colormap(ha(i),Color_KGE)
    clim([-0.4,1])
    c = colorbar;
    set(c,'ticks',(-0.2:0.2:0.8),'ticklength',0.015,'fontsize',13,'FontName','Times New Roman')
    h_title = title(cell2mat(title_name(i)),'FontSize',16,'FontName','Time New Roman',...
        'HorizontalAlignment','left','Units', 'normalized');
    set(h_title, 'Position',[0.005, 1.025, 0]);
    txt1 = ['Median: ',num2str(round(median(kge(:,i),'omitnan'),2));];
    text(0.005, 0.05, txt1, 'Units', 'normalized', 'fontsize', 15, 'FontName', 'Time New Roman', 'Color', 'r','FontWeight','bold');
    for j = 1:size(Basin_SHP,1)
        m_line(Basin_SHP(j).Lon,Basin_SHP(j).Lat,'color',[0.5 0.5 0.5, 0.5], 'linewidth', 0.6);
    end
end
print(figure(1),'./Fig3_Streamflow_KGE','-dpng','-r300')

% Calculate summary statistics (Median) for KGE and RMSE across all stations
mid_kge = median(kge,1,'omitnan');
mid_rmse = median(rmse,1,'omitnan');

fprintf('%-15s | %-10s | %-10s | %-15s | %-15s\n', ...
    'Model Name', 'Median KGE', 'Median RMSE', 'KGE Increase', 'RMSE Improve (%)');
model_names = {'CSSPv3','CSSPv2','DBH','H08','LPJML','PCR-GLOBWB','WATERGAP','CLM4'};
for i = 1:8
    if i == 1
        kge_inc_str = '-';
        rmse_imp_str = '-';
    else
        kge_inc = mid_kge(1) - mid_kge(i);
        kge_inc_str = num2str(round(kge_inc, 2));
        rmse_imp = (mid_rmse(i) - mid_rmse(1)) / mid_rmse(i) * 100;
        rmse_imp_str = [num2str(round(rmse_imp, 1)), '%'];
    end
    fprintf('%-15s | %-10.2f | %-10.2f | %-15s | %-15s\n', ...
        model_names{i}, mid_kge(i), mid_rmse(i), kge_inc_str, rmse_imp_str);
end
