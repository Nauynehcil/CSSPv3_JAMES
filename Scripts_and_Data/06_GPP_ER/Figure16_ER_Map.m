clc,clear
load ../mask.mat
load ./FluxSat.mat
load ('./CSSP_output.mat','CSSP_MEAN');
load ('./CSSP_output.mat','CSSP_LAT');
load ('./MsTMIP_output.mat','MsTMIP_MEAN');
load ('./MsTMIP_output.mat','MsTMIP_LAT');

M = [243,243,243; 252,248,187; 254,191,136; 246,137,109; 226,085,105; 184,056,116; 131,037,127; 085,020,118; 044,016,076]/255;
Heat_Color = [255,255,255; 235,235,235; 252,248,187; 254,191,136; 246,137,109; 226,085,105; 184,056,116; 131,037,127; 085,020,118; 044,016,076]/255;

fig_map(:,:,1) = CSSP_MEAN.ER_mon(:,:,1);
fig_map(:,:,2) = median(MsTMIP_MEAN.ER_mon,4,'omitnan');

CSSPv3_LAT(:,1) = CSSP_LAT.ER_mon(:,1);
MsTMIP_mid_LAT(:,1) = prctile(MsTMIP_LAT.ER_mon,50,2);
MsTMIP_upper_LAT(:,1) = prctile(MsTMIP_LAT.ER_mon,75,2);
MsTMIP_lower_LAT(:,1) = prctile(MsTMIP_LAT.ER_mon,25,2);

figure(1)
set(gcf,'Units','centimeters','Position',[5,5,32.5,20.5],'color','white');
ha = tight_subplot(2,1,[.025 .05],[.025 .025],[.055 .325]);
ma = tight_subplot(1,1,[.0 .05],[.075 .065],[.725 .045]);
title_name = {'(a) CSSPv3 ER','(b) MsTMIP Median ER','(c) Mean Latitudinal Profile of ER'};

for i = 1:2
    axes(ha(i))
    m_proj('Equidistant Cylindrical','long',[-180 180],'lat',[-60 90]);hold on
    m_pcolor(xlon,xlat,fig_map(:,:,i));hold on
    m_grid('xtick',(-180:30:180),'ytick',(-60:30:90),'tickdir','out','FontSize',12,'FontName','Times New Roman');

    clim([0 270])
    colormap(ha(i),M)
    c=colorbar;
    set(c,'ticks',(0:30:240),'ticklength',0.015,'fontsize',12,'FontName','Times New Roman')
    text(1, 1.025,'Unit: gC/m^2/month', 'Units', 'normalized', ...
        'FontSize', 14, 'FontName', 'Times New Roman', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

    h_title = title(cell2mat(title_name(i)),'FontSize',14,'FontName','Time New Roman',...
        'HorizontalAlignment','left','Units', 'normalized');
    set(h_title, 'Position',[0.005, 1.025, 0]);
    for j = 1:size(Basin_SHP,1)
        m_line(Basin_SHP(j).Lon,Basin_SHP(j).Lat,'color',[0.5 0.5 0.5], 'linewidth', 1);
    end
end

axes(ma(1))
LAT(:,1) = xlat(1,:);
valid_idx = ~isnan(LAT) & ~isnan(MsTMIP_lower_LAT) & ~isnan(MsTMIP_upper_LAT) & ~isnan(MsTMIP_mid_LAT) ;
LAT_valid = LAT(valid_idx);
CSSPv3_valid = CSSPv3_LAT(valid_idx,1);
MsTMIP_mid_valid = MsTMIP_mid_LAT(valid_idx);
MsTMIP_lower_valid = MsTMIP_lower_LAT(valid_idx);
MsTMIP_upper_valid = MsTMIP_upper_LAT(valid_idx);

IQR = fill([MsTMIP_lower_valid; flipud(MsTMIP_upper_valid)], ...
    [LAT_valid; flipud(LAT_valid)], ...
    [186,215,231]/255,'FaceAlpha',0.5,'EdgeColor', 'none',...
    'DisplayName', 'ALL IQR (25%-75%)');
hold on
fig_lat(1) = plot(CSSPv3_valid, LAT_valid, 'LineWidth', 1.5, 'Color', [246,001,001]/255);
fig_lat(2) = plot(MsTMIP_mid_valid, LAT_valid, 'LineWidth', 1.5, 'Color', [000,034,164]/255);

hold off;
set(gca,'XLim',[0 280], 'XTick', 0:80:280, 'fontsize',13, 'FontName','Times New Roman');
set(gca,'YLim',[-40 78]);
Y = -90:15:90;
set(gca,'YTick',Y, 'fontsize',12, 'FontName','Times New Roman');
lat_str = arrayfun(@(y) sprintf('%d°%s', abs(y), char(78*(y>0) + 83*(y<0))), Y, 'UniformOutput', false);
set(gca,'YTickLabel',lat_str, 'FontSize',13, 'Fontname','Times New Roman');
xlabel('[g.C/m^2/month]','FontSize',13,'Fontname', 'Time New Roman')
m_title = title(cell2mat(title_name(3)),'FontSize',14,'FontName','Time New Roman',...
    'HorizontalAlignment','left', 'Units', 'normalized');
set(m_title, 'Position', [0.005, 1.005, 0]); 
l = legend([fig_lat(1), fig_lat(2), IQR], ...
    {'CSSPv3','Ensemble Median', 'Ensemble IQR'}, ...
    'FontSize', 13, 'FontName', 'Times New Roman');
box off
ax = axes('Position',get(gca,'Position'),'XAxisLocation','top',...
    'YAxisLocation','right','Color','none','XColor','k','YColor','k'); 
set(ax,'XTick', [],'YTick', []); 
print(figure(1),'./Fig3_ER','-dpng','-r300')
