clc,clear
load ../mask.mat
load ./Metric_Data.mat

% Custom colormap for KGE: transitioning from Blue (Low/Poor) to Dark Red (High/Good)
Color_KGE = [069,117,180; 224,243,248; 255,245,229; 253,209,155; 252 140 089; 216,051,033; 135,000,000]/255;

% Extract KGE for CSSPv3
CSSP_kge = kge(:,:,1);

% title name
title_name = {'(a) KGE of CSSPv3-Simulated SM at Surface Layer',...
              '(b) KGE of CSSPv3-Simulated SM at Middle Layer',...
              '(c) KGE of CSSPv3-Simulated SM at Deep Layer'};

% Binning Parameters
res = 1; % Grid resolution (degrees)
lon_edges = -180:res:180;
lat_edges = -56:res:84;

figure(1)
set(gcf,'Units','centimeters','Position',[5,1.5,25.5,32.5],'color','white')
ha = tight_subplot(3,1,[.05 .05],[.05 .05],[.1 .1]);
for i = 1:3
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
    binned_kge = accumarray([binY(valid), binX(valid)], CSSP_kge(valid, i), sz, @(x) median(x, 'omitnan'), NaN);
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
    txt1 = ['Median: ',num2str(round(median(CSSP_kge(:,i),'omitnan'),2));];
    text(0.005, 0.05, txt1, 'Units', 'normalized', 'fontsize', 15, 'FontName', 'Time New Roman', 'Color', 'r','FontWeight','bold');
    text(1.0, 1.025,[num2str(count(i)),' Merged Sites'], 'Units', 'normalized', ...
        'FontSize', 16, 'FontName', 'Times New Roman', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

    for j = 1:size(Basin_SHP,1)
        m_line(Basin_SHP(j).Lon,Basin_SHP(j).Lat,'color',[0.5 0.5 0.5, 0.5], 'linewidth', 0.6);
    end
end
print(figure(1),'./Fig7_SM_CSSP_KGE','-dpng','-r300')





