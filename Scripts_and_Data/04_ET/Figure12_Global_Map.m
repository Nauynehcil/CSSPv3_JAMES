clc,clear
load ./CSSP.mat
load ./Models.mat

ALL_LAT_mon = cat(1, GLEAM_LAT_mon, ERA5_LAT_mon, ERA5_LAND_LAT_mon, ISIMIPA_LAT_mon', MsTMIP_LAT_mon');
ALL_MEAN_mon = cat(3, MEAN_GLEAM_mon, MEAN_ERA5_mon, MEAN_ERA5_LAND_mon, MEAN_ISIMIPA_mon, MEAN_MsTMIP_mon);
clearvars -except ALL* *CSSP*

load ../mask.mat
M = [243,243,243; 252,248,187; 254,191,136; 246,137,109; 226,085,105; 184,056,116; 131,037,127; 085,020,118; 044,016,076]/255;
N = [027,059,112; 036,107,176; 074,143,190; 123,178,212; 186,215,231; 221,232,238; 242,242,242;...
    242,242,242; 246,228,218; 248,204,183; 240,158,123; 216,099,079; 185,033,045; 112,012,034]/255;

fig_map(:,:,1) = MEAN_CSSPv3_mon;
fig_map(:,:,2) = median(ALL_MEAN_mon,3,'omitnan');
ALL_mid_LAT(:,1) = prctile(ALL_LAT_mon,50,1);
ALL_upper_LAT(:,1) = prctile(ALL_LAT_mon,75,1);
ALL_lower_LAT(:,1) = prctile(ALL_LAT_mon,25,1);
title_name = {'(a) CSSPv3 ET','(b) Ensemble Median ET','(c) Mean Latitudinal Profile of ET'};

figure(1)
set(gcf,'Units','centimeters','Position',[5,5,32.5,20.5],'color','white');
ha = tight_subplot(2,1,[.025 .05],[.025 .025],[.055 .325]);
ma = tight_subplot(1,1,[.0 .05],[.075 .065],[.725 .045]);

for i = 1:2
    axes(ha(i))
    m_proj('Equidistant Cylindrical','long',[-180 180],'lat',[-60 90]);hold on
    m_pcolor(xlon,xlat,fig_map(:,:,i));hold on
    m_grid('xtick',(-180:30:180),'ytick',(-60:30:90),'tickdir','out','FontSize',12,'FontName','Times New Roman');
    clim([0 180])
    colormap(ha(i),M)
    c=colorbar;
    set(c,'ticks',(0:20:160),'ticklength',0.015,'fontsize',12,'FontName','Times New Roman')
    text(1.0, 1.025,'Unit: mm/month', 'Units', 'normalized', ...
        'FontSize', 15, 'FontName', 'Times New Roman', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

    h_title = title(cell2mat(title_name(i)),'FontSize',15,'FontName','Time New Roman',...
        'HorizontalAlignment','left','Units', 'normalized');
    set(h_title, 'Position',[0.005, 1.025, 0]);
    for j = 1:size(Basin_SHP,1)
        m_line(Basin_SHP(j).Lon,Basin_SHP(j).Lat,'color',[0.5 0.5 0.5], 'linewidth', 1);
    end
end

axes(ma(1))
LAT(:,1) = xlat(1,:);
valid_idx = ~isnan(LAT) & ~isnan(ALL_lower_LAT) & ~isnan(ALL_upper_LAT) & ~isnan(ALL_mid_LAT) ;
LAT_valid = LAT(valid_idx);
CSSPv3_valid = CSSPv3_LAT_mon(valid_idx);
CSSPv2_valid = CSSPv2_LAT_mon(valid_idx);
ALL_mid_valid = ALL_mid_LAT(valid_idx);
ALL_lower_valid = ALL_lower_LAT(valid_idx);
ALL_upper_valid = ALL_upper_LAT(valid_idx);
IQR = fill([ALL_lower_valid; flipud(ALL_upper_valid)], ...
    [LAT_valid; flipud(LAT_valid)], ...
    [186,215,231]/255,'FaceAlpha',0.5,'EdgeColor', 'none',...
    'DisplayName', 'ALL IQR (25%-75%)');
hold on
fig_lat(2) = plot(CSSPv2_valid, LAT_valid, 'LineWidth', 1.5, 'Color', [077,175,074]/255);
fig_lat(1) = plot(CSSPv3_valid, LAT_valid, 'LineWidth', 1.5, 'Color', [246,001,001]/255);
fig_lat(3) = plot(ALL_mid_valid, LAT_valid, 'LineWidth', 1.5, 'Color', [000,034,164]/255);
hold off;
set(gca,'XLim',[0 130], 'XTick', 0:50:150, 'fontsize',13, 'FontName','Times New Roman');
set(gca,'YLim',[-40 78]);
Y = -90:15:90;
set(gca,'YTick',Y, 'fontsize',12, 'FontName','Times New Roman');
lat_str = arrayfun(@(y) sprintf('%d°%s', abs(y), char(78*(y>0) + 83*(y<0))), Y, 'UniformOutput', false);
set(gca,'YTickLabel',lat_str, 'FontSize',13, 'Fontname','Times New Roman');
xlabel('[mm/month]','FontSize',14,'Fontname', 'Time New Roman')
m_title = title(cell2mat(title_name(3)),'FontSize',15,'FontName','Time New Roman',...
    'HorizontalAlignment','left', 'Units', 'normalized');
set(m_title, 'Position', [0.005, 1.005, 0]); 
l = legend([fig_lat(1), fig_lat(2), fig_lat(3), IQR], ...
    {'CSSPv3', 'CSSPv2', 'Ensemble Median', 'Ensemble IQR'}, ...
    'FontSize', 13, 'FontName', 'Times New Roman');
box off
ax = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
    'YAxisLocation','right','Color','none','XColor','k','YColor','k');
set(ax,'XTick', [],'YTick', []);
print(figure(1),'./Fig12_FEVPA_MAP','-dpng','-r300')
