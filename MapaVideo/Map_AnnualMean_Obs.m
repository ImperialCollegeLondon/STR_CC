% This file is to map averaged annual precipitation from
% several observation datasets
% Including:
%        GEAR
%        HAD-UK GRID
%        NIMROD RADAR
%
% @ Yuting Chen
% Imperial College London
% yuting.chen17@imperial.ac.uk

clear;clc

REGIONS = REGIONS_info();
region = REGIONS.UK;% Scotland;% Westuk;% wales

%% hadukgrid
[HadUK,E_had,N_had,HadUKyr] = getHadUK('tas');

figure;
setFigureProperty('Paper');
mapIt(E_had,N_had,HadUKyr)
text(400,1100,sprintf('HadUKGrid(1981-2000)'),'horizontalalignment','center')




%% hadukgrid
[HadUK,E_had,N_had,HadUKyr] = getHadUK('pr');

figure;
setFigureProperty('Paper');
mapIt(E_had,N_had,HadUKyr)
text(400,1100,sprintf('HadUKGrid(1981-2000)'),'horizontalalignment','center')

filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP';
filename = [filePath,filesep,'HadUKGrid_annual'];
savePlot(filename,'units','centimeters','XYWH',[5,0,8,11],'needreply','Y','onlyPng',false);

close all

%% gear
[GEAR,X_coor,Y_coor,DIST] = getAverageGEAR(region,'month');
gear1km = squeeze(nansum(GEAR,1));
gear1km(gear1km==0) = NaN;

figure;
setFigureProperty('Paper');
mapIt(X_coor,Y_coor,gear1km);
text(400,1100,sprintf('GEAR(1981-2000)'),'horizontalalignment','center')


filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP';
filename = [filePath,filesep,'GEAR1980_1999'];
savePlot(filename,'units','centimeters','XYWH',[5,0,8,11],'needreply','Y','onlyPng',false);

%% nimrod radar
[RAD,E_rad,N_rad] = getMonthRadar(region);
in = getTrimTag('unit','km','product','radar2.2');
rad = squeeze(nansum(eomday(1999,1:12)'.*RAD,1))*24;
rad(~in) = NaN;

figure;
setFigureProperty('Paper');
mapIt(E_rad,N_rad,rad);
text(400,1100,sprintf('Radar(2007-2018)'),'horizontalalignment','center')

filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP';
filename = [filePath,filesep,'Radar2007_2018'];
savePlot(filename,'units','centimeters','XYWH',[5,0,8,11],'needreply','Y','onlyPng',false);

function mapIt(E,N,Z)
UKMap = getUKMap();
pcolor(E,N,Z);shading flat
cptcmap('precip_annualMeanUK', 'mapping','direct');%,'ncol',20);
axis off
c = colorbar('location','Manual', 'position', [0.65 0.45 0.02 0.3],'fontsize',12);
c.Ruler.TickLabelFormat='%gmm';
hold on; plot(UKMap.borderE/1000,UKMap.borderN/1000,'-','color',[0.5,0.5,0.5],'linewidth',1);
ylim([0,1200]);caxis([0,3100]);xlim([0,800]);
axis equal
xlim([0,800]);
end


