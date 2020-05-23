% --------------------------------------------------------- %
% this file is to plot all regions (110Km*110Km per each) of interest over UK as rectangles.
% and then shown in UK map with a background of GEAR rainfall data
% --------------------------------------------------------- %


clear;clc
close all

REGIONS = REGIONS_info();
region = REGIONS.UK;

% Get Result from RADAR-Month
[RAIN,X_coor,Y_coor] = getMonthRadar(region);
in = getTrimTag('unit','km','product','radar2.2');
RAIN = squeeze(nansum(eomday(1999,1:12)'.*RAIN,1))*24;
RAIN(~in) = NaN;
rain2d2km = RAIN;

%%
UKMap = load('H:\CODE_MATLAB\SpatialTemporalDATA\shapeFileFolder\UKBorderGrid.mat');

figure
setFigureProperty('Paper');
hold on

pcolor(X_coor,Y_coor,rain2d2km);shading flat

cptcmap('precip_annualMeanUK', 'mapping','direct');%,'ncol',20);
axis off
% c = colorbar;%('location','Manual', 'position', [0.1 0.1 0.02 0.81],'fontsize',15);
% c = colorbar('southoutside');
c = colorbar('location','Manual', 'position', [0.65 0.45 0.02 0.3],'fontsize',12);

c.Ruler.TickLabelFormat='%gmm';
c.FontSize = 12;
plot(UKMap.borderE/1000,UKMap.borderN/1000,'-','linewidth',1,'color',[0.5,0.5,0.5]);
ylim([0,1200]);caxis([0,3100]);xlim([0,800]);

% get several regions of interest
REGIONS = REGIONS_info();

plotOneRe = @(region)rectangle('position',[region.minE,region.minN,110,110],...
    'linewidth',1,'Linestyle','-');
tagOneRe = @(region)text(region.minE+55,region.minN+55,region.Name,...
    'fontsize',8,'horizontalalignment','center','fontweight','bold',...
    'background',[0.9 0.9 0.9 0.3]);

plotOneRe(REGIONS.SWestuk);
plotOneRe(REGIONS.Westuk);
plotOneRe(REGIONS.Scotland);
plotOneRe(REGIONS.London);

tagOneRe(REGIONS.SWestuk);
tagOneRe(REGIONS.Westuk);
tagOneRe(REGIONS.Scotland);
tagOneRe(REGIONS.London);

hold off
axis equal
xlim([0,800]);

fileName = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\regionsLocationMap';
savePlot(fileName,'units','centimeters','XYWH',[5,0,8,11],'needreply','Y');




