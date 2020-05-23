% ----------------------------------------------------------------------- %
% THIS FILE IS TO INVESTIGATE IF CLIMATE MODELS PRODUCE WELL ON AVERAGE
% RAINFALL VALUE
% 
% DATA USED:
%           Climate Data: UKCP18 CPM 2.2
%           RainGauge Data: Monthly GEAR from CEH
%           Radar Data: MetOffice NIMROD Radar Composite 5min
%           Topo Data: Ordance Survey TERRAIN50 DTM
% 
% @ Yuting CHEN
% yuting.chen17@imperial.ac.uk
% Update: 2020.01.27
% ----------------------------------------------------------------------- %



clear;clc
setFigureProperty('Paper');
close all

REGIONS = REGIONS_info();
region = REGIONS.UK;%Scotland;% Westuk;% wales

%% Get All data required: GEAR-Month, HadUK-Grid, RADAR, CPM2.2

% Result from GEAR-Month
[GEAR,X_coor,Y_coor,DIST] = getAverageGEAR(region,'month');


% UKMap = load('H:\CODE_MATLAB\SpatialTemporalDATA\shapeFileFolder\UKBorderGrid.mat');
% in = inpolygon(X_coor,Y_coor,UKMap.borderE/1000,UKMap.borderN/1000);
% GEAR(:,~in) = NaN; 

% [gearIJ] = find(DIST>1000);
% GEAR(:,gearIJ) = NaN;


% Result from HadUK-Grid
[HadUK,E_had,N_had,HadUKyr] = getHadUK();

% Result from RADAR
[RAD,E_rad,N_rad] = getMonthRadar(region);

% Result from CPM
[CPM,E,N] = getMonthCPM(region);% [mon,ensNo]


%% PLOT MEAN SEASONAL CYCLE OF PRECIPITATION (IPDD 2018 Evaluation of climate models Fig 9.38)

ha = tight_subplot(2,2,[.1 .1],[.10 .10],[.10 .30]);

set(gcf,'units','points','position',[150,0,700,500]);


axes(ha(1))
plotSeasonality(REGIONS.London,CPM,E,N,GEAR,X_coor,Y_coor,DIST,RAD,E_rad,N_rad)
axes(ha(2))
plotSeasonality(REGIONS.SWestuk,CPM,E,N,GEAR,X_coor,Y_coor,DIST,RAD,E_rad,N_rad)
axes(ha(3))
plotSeasonality(REGIONS.Westuk,CPM,E,N,GEAR,X_coor,Y_coor,DIST,RAD,E_rad,N_rad)
axes(ha(4))
plotSeasonality(REGIONS.Scotland,CPM,E,N,GEAR,X_coor,Y_coor,DIST,RAD,E_rad,N_rad)

legend boxon

filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP';
filename = [filePath,filesep,'regional_monthly'];

savePlot(filename,'XYWH',[150,0,700,500],'needreply','Y');


%% PLOT WHOLE PATTERN (Annual)
close all
gear1km = squeeze(nansum(GEAR,1));
gear1km(gear1km==0) = NaN;

in = getTrimTag('unit','km','product','radar2.2');
rad = squeeze(nansum(eomday(1999,1:12)'.*RAD,1))*24;
rad(~in) = NaN;

haduk = HadUKyr;

in = getTrimTag('unit','km','product','cpm2.2');
cpm12 = [];
for ensNo = 1:12
    cpm0 = nansum(cat(3,CPM{:,ensNo}).*reshape(eomday(1999,1:12),[1,1,12]), 3)*24;
    cpm0(~in) = NaN;
    cpm12 = cat(3,cpm12,reshape(cpm0,[size(cpm0),1]));
end
cpm = nanmedian(cpm12,3);


f_1to2d2 = @(x)imresize(imresize(x,5,'nearest'),1/11,'box');
g2 = f_1to2d2(gear1km);
x2 = f_1to2d2(X_coor);
y2 = f_1to2d2(Y_coor);
Gear2d2 = griddata(x2,y2,g2,E,N,'nearest');
bias_gear = (cpm-Gear2d2)./Gear2d2;
bias_gear_sig = testSig(permute(cpm12-Gear2d2,[3,1,2]));
bias_gear(bias_gear_sig==0) = 0;%
g2 = f_1to2d2(haduk);
x2 = f_1to2d2(E_had);
y2 = f_1to2d2(N_had);
had2d2 = griddata(x2,y2,g2,E,N,'nearest');
bias_had = (cpm-had2d2)./had2d2;
bias_had_sig = testSig(permute(cpm12-had2d2,[3,1,2]));
bias_had(bias_had_sig==0) = 0;%

plotAllAvera([],X_coor,Y_coor,gear1km,E_rad,N_rad,rad,...
    E_had,N_had,HadUKyr,E,N,cpm,bias_gear,bias_had)

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

filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\bias';
filename = [filePath,filesep,'gear_radar_cpm_bias'];
savePlot(filename,'units','centimeters','XYWH',[5,0,24,22],'needreply','Y');


%% PLOT JJA pattern
close all
mons = 6:8;
gear1km = squeeze(nansum(GEAR(mons,:,:),1));
gear1km(gear1km==0) = NaN;

in = getTrimTag('unit','km','product','radar2.2');
rad = squeeze(nansum(eomday(1999,mons)'.*RAD(mons,:,:),1))*24;
rad(~in) = NaN;

haduk = squeeze(nansum(HadUK(:,:,mons),3));
haduk(haduk==0) = NaN;

in = getTrimTag('unit','km','product','cpm2.2');
cpm12 = [];
for ensNo = 1:12
    cpm0 = nansum(cat(3,CPM{mons,ensNo}).*reshape(eomday(1999,mons),...
        [1,1,length(mons)]), 3)*24;
    cpm0(~in) = NaN;
    cpm12 = cat(3,cpm12,reshape(cpm0,[size(cpm0),1]));
end
cpm = nanmedian(cpm12,3);

f_1to2d2 = @(x)imresize(imresize(x,5,'nearest'),1/11,'box');
g2 = f_1to2d2(haduk);
x2 = f_1to2d2(E_had);
y2 = f_1to2d2(N_had);
had2d2 = griddata(x2,y2,g2,E,N,'nearest');
bias_had = (cpm-had2d2)./had2d2;
bias_had_sig = testSig(permute(cpm12-had2d2,[3,1,2]));
bias_had(bias_had_sig==0) = 0;%

g2 = f_1to2d2(gear1km);
x2 = f_1to2d2(X_coor);
y2 = f_1to2d2(Y_coor);
Gear2d2 = griddata(x2,y2,g2,E,N,'nearest');
bias_gear = (cpm-Gear2d2)./Gear2d2;
bias_gear_sig = testSig(permute(cpm12-Gear2d2,[3,1,2]));
bias_gear(bias_gear_sig==0) = 0;%

% see radar;
% f_1to2d2 = @(x)imresize(imresize(x,5,'nearest'),1/11,'box');
% g2 = rad;
% x2 = E_rad;
% y2 = N_rad;
% rad2d2 = griddata(x2,y2,g2,E,N,'nearest');
% bias = (cpm-rad2d2)./rad2d2;

sf = 365/sum(eomday(1999,mons));

plotAllAvera([],X_coor,Y_coor,gear1km*sf,E_rad,N_rad,rad*sf,...
    E_had,N_had,HadUKyr*sf,E,N,cpm*sf,bias_gear,bias_had)


filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\bias';
filename = [filePath,filesep,'gear_radar_cpm_bias_JJA'];
savePlot(filename,'units','centimeters','XYWH',[5,0,24,22],'needreply','Y');
close all

%% PLOT THESE THREE AND BIAS FOR EACH MONTH SEPERATELY

for mon = 1:12
    
    gear1km = squeeze(GEAR(mon,:,:))/eomday(1990,mon)*365;
    
    in = getTrimTag('unit','km','product','radar2.2');
    rad = squeeze(RAD(mon,:,:))*24*365;
    rad(~in) = NaN;
    
    haduk = squeeze(HadUK(:,:,mon))/eomday(1990,mon)*365;
    haduk(haduk==0) = NaN;
    
    in = getTrimTag('unit','km','product','cpm2.2');
    cpm = nanmean(cat(3,CPM{mon,:}), 3)*24*365;
    cpm(~in) = NaN;
    
    
    g2 = f_1to2d2(gear1km);
    x2 = f_1to2d2(X_coor);
    y2 = f_1to2d2(Y_coor);
    Gear2d2 = griddata(x2,y2,g2,E,N,'nearest');
    bias_gear = (cpm-Gear2d2)./Gear2d2;
    
    g2 = f_1to2d2(haduk);
    x2 = f_1to2d2(E_had);
    y2 = f_1to2d2(N_had);
    had2d2 = griddata(x2,y2,g2,E,N,'nearest');
    bias_had = (cpm-had2d2)./had2d2;

    %%
    plotAllAvera(mon,X_coor,Y_coor,gear1km,E_rad,N_rad,rad,...
        E_had,N_had,haduk,E,N,cpm,bias_gear,bias_had)
    
    filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\bias';
    filename = [filePath,filesep,'gear_radar_cpm_bias_Mon',sprintf('%02d',mon)];
    savePlot(filename,'units','centimeters','XYWH',[5,0,24,22],'needreply','N');
    close all
    
end

%% AUXILLARY FUNCTION
function [isSigChange] = testSig(d)
% d: 3d# [series, dim1, dim2]
% 1: different 
% 0: not different.

isSigChange = squeeze(ttest(d,0,'alpha',0.1));

end

function plotSeasonality(region,CPM,E,N,GEAR,X_coor,Y_coor,DIST,RAD,E_rad,N_rad)

%CPM
for ensNo = 1:12;
    [i1,j1] = getRegionIJ(E,N,region.minE,region.minN);
    [i2,j2] = getRegionIJ(E,N,region.minE+region.dx*(region.dimE-1),...
        region.minN+region.dx*(region.dimN-1));
    CPMthis = cat(3,CPM{:,ensNo});CPMthis = CPMthis(i1:i2,j1:j2,:);
    CPM_mon(ensNo,:) = 24*squeeze(nanmean(reshape(CPMthis,[],12),1));
end
hsimran = fill([1:12,12:-1:1],[nanmin(CPM_mon,[],1),flip(nanmax(CPM_mon,[],1))],'r',...
    'LineStyle','none');hold on;
alpha(0.3)
hsim = plot(1:12,nanmean(CPM_mon,1),'--','color','r');hold on


% GEAR
GEAR = permute(GEAR,[1,3,2]);
X_coor = permute(X_coor,[2,1]);
Y_coor = permute(Y_coor,[2,1]);
DIST = permute(DIST,[2,1]);
GEAR(GEAR == 0) = NaN;
[i1,j1] = getRegionIJ(X_coor,Y_coor,region.minE,region.minN);
[i2,j2] = getRegionIJ(X_coor,Y_coor,region.minE+region.dx*(region.dimE-1),...
    region.minN+region.dx*(region.dimN-1));
GEARthis = GEAR(:,i1:i2,j1:j2);
GEARthis = permute(GEARthis,[2,3,1]);

GEAR_mon = squeeze(nanmean(reshape(GEARthis,[],12),1))./eomday(1990,1:12);
hobs_gear = plot(1:12,GEAR_mon,'-','color','k','linewidth',3);hold on

% RAD
in = getTrimTag('unit','km','product','radar2.2');
in = repmat(reshape(in,[1,size(in)]),[12,1,1]);
RAD(~in) = NaN;
[i1,j1] = getRegionIJ(E_rad,N_rad,region.minE,region.minN);
[i2,j2] = getRegionIJ(E_rad,N_rad,region.minE+region.dx*(region.dimE-1),...
    region.minN+region.dx*(region.dimN-1));
RADthis = RAD(:,i1:i2,j1:j2);
RADthis = permute(RADthis,[2,3,1]);

RAD_mon = 24*squeeze(nanmean(reshape(RADthis,[],12),1));
hobs_rad = plot(1:12,RAD_mon,'--','color','k','linewidth',3);

xlim([1,12]);
xticks([1:12])
xticklabels(getMonthName(1:12,1,1))
YLIM = ylim;
text(1.2,YLIM(2),region.Name,'verticalAlignment','top','horizontalAlignment','left')


legend([hsimran,hsim,hobs_gear,hobs_rad],{'EnsembleRange','CPM2.2 mean','GEAR','RAD'})
legend off
box on;
set(gca,'linewidth',2)
end



function plotAllAvera(mon,X_coor,Y_coor,gear1km,E_rad,N_rad,rad,...
    E_had,N_had,HadUK,E,N,cpm,bias_gear,bias_had)

UKMap = load('H:\CODE_MATLAB\SpatialTemporalDATA\shapeFileFolder\UKBorderGrid.mat');

setFigureProperty('Meeting');
ha = tight_subplot(2,3,[.05 .05],[.05 .05],[.15 .15]);


axes(ha(3))
Z = HadUK;
pcolor(E_had,N_had,Z);shading flat
cptcmap('precip_annualMeanUK', 'mapping','direct');%,'ncol',20);
xlim([0,700]);ylim([0,1200]);caxis([0,3100])
axis off
text(200,1100,sprintf('HadUKGrid(1981-2000)%s',getMonthName(mon)))
c = colorbar('location','Manual', 'position', [0.1 0.1 0.02 0.81],'fontsize',15);
% c = colorbar('southoutside');
c.Ruler.TickLabelFormat='%gmm';
c.FontSize = 12;
hold on; plot(UKMap.borderE/1000,UKMap.borderN/1000,'-','color',[0.5,0.5,0.5]);


axes(ha(1))
pcolor(X_coor,Y_coor,gear1km);shading flat
cptcmap('precip_annualMeanUK', 'mapping','direct');%,'ncol',20);
xlim([0,700]);ylim([0,1200]);caxis([0,3100])
axis off
text(200,1100,sprintf('GEAR %s',getMonthName(mon)))

hold on; plot(UKMap.borderE/1000,UKMap.borderN/1000,'-','color',[0.5,0.5,0.5]);

axes(ha(2))
Z = rad;
pcolor(E_rad,N_rad,Z);shading flat
cptcmap('precip_annualMeanUK', 'mapping','direct');%,'ncol',20);
xlim([0,700]);ylim([0,1200]);caxis([0,3100])
axis off
text(200,1100,sprintf('RAD(2007-2018) %s',getMonthName(mon)))

hold on; plot(UKMap.borderE/1000,UKMap.borderN/1000,'-','color',[0.5,0.5,0.5]);

axes(ha(4))
Z = cpm;
pcolor(E,N,Z);shading flat
cptcmap('precip_annualMeanUK', 'mapping','direct');%,'ncol',20);
xlim([0,700]);ylim([0,1200]);caxis([0,3100])
axis off
text(200,1100,sprintf('CPM(1980-2000) %s',getMonthName(mon)))
c = colorbar('location','Manual', 'position', [0.1 0.1 0.02 0.81],'fontsize',15);
% c = colorbar('southoutside');
c.Ruler.TickLabelFormat='%gmm';
c.FontSize = 12;
hold on; plot(UKMap.borderE/1000,UKMap.borderN/1000,'-','color',[0.5,0.5,0.5]);


axes(ha(5))
pcolor(E,N,bias_gear*100);shading flat
cptcmap('diff_darkBlue_darkRed', 'mapping','direct','ncol',11);
caxis([-110,110])
hold on; plot(UKMap.borderE/1000,UKMap.borderN/1000,'k-');
xlim([0,700])
ylim([0,1200])
axis off
text(200,1100,sprintf('Bias (CPM-GEAR)%s',getMonthName(mon)))

axes(ha(6))
pcolor(E,N,bias_had*100);shading flat
cptcmap('diff_darkBlue_darkRed', 'mapping','direct','ncol',11);%diff_4_bias
caxis([-110,110])
hold on; plot(UKMap.borderE/1000,UKMap.borderN/1000,'k-');
xlim([0,700])
ylim([0,1200])
axis off
text(200,1100,sprintf('Bias (CPM-HadUK)%s',getMonthName(mon)))

c = colorbar('location','Manual', 'position', [0.85 0.1 0.02 0.30]);
c.Ticks = [-100,-30,0,30,100];
c.Ruler.TickLabelFormat='%g%%';
c.FontSize = 12;



end




