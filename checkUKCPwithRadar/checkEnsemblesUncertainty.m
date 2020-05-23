%
% This file is to check CPM ensembles.
% mainly focused the uncertainty in these 12 ensemble simulation
%
% Including:
%
% Interannual variability
%
% @ Yuting Chen
% yuting.chen17@imperial.ac.uk
% Imperial College London
%

clear;clc

REGIONS = REGIONS_info();
region = REGIONS.UK;%Scotland;% Westuk;% wales
ENSEMBLENO=getEnsNos();


%% scheme 1: annual change against year 
% <1980-1999>



%% get CPM + GEAR (whole UK)
[CPM,E,N] = getYearCPM(region);% [mon,ensNo]
CPM = permute(CPM,[1,4,2,3]);%[mon,yr,E,N]
in = getTrimTag('unit','km','product','cpm2.2');
in = reshape(in,[1,1,size(CPM,3:4)]);
CPM = CPM.*in;
CPM_size = size(CPM);
CPM(CPM==0) = NaN;


GEAR = struct;
[GEAR.val,GEAR.X_coor,GEAR.Y_coor,GEAR.dist] = getAverageGEAR(region,'year');
GEAR.val = permute(GEAR.val,[1,3,2]);
GEAR.X_coor = permute(GEAR.X_coor,[2,1]);
GEAR.Y_coor = permute(GEAR.Y_coor,[2,1]);
GEAR.dist = permute(GEAR.dist,[2,1]);
GEAR.val(GEAR.val == 0) = NaN;

%% Get a regional result
figure;
ha = tight_subplot(3,2,[.1 .1],[.15 .15],[.10 .10]);

axes(ha(1))
[sim,obs,STATS,GEARthis,CPMthis]=getRegionPlot(REGIONS.London,E,N,CPM,GEAR);

axes(ha(2))
[sim,obs,STATS(2,:),GEARthis,CPMthis]=getRegionPlot(REGIONS.SWestuk,E,N,CPM,GEAR);

axes(ha(3))
[sim,obs,STATS(3,:),GEARthis,CPMthis]=getRegionPlot(REGIONS.Westuk,E,N,CPM,GEAR);

axes(ha(4))
[sim,obs,STATS(4,:),GEARthis,CPMthis]=getRegionPlot(REGIONS.Scotland,E,N,CPM,GEAR);

subplot(3,2,[5,6])
[sim,obs,STATS(6,:),GEARthis,CPMthis]=getRegionPlot(REGIONS.smallUK,E,N,CPM,GEAR);



%%

plot_ENESEMBLE_MONTHLY_MEAN();

function [sim,obs,STATS,GEARthis,CPMthis] = getRegionPlot(region,E,N,CPM,GEAR)

[i1,j1] = getRegionIJ(E,N,region.minE,region.minN);
[i2,j2] = getRegionIJ(E,N,region.minE+region.dx*(region.dimE-1),...
    region.minN+region.dx*(region.dimN-1));
CPMthis = CPM(:,:,i1:i2,j1:j2);

[i1,j1] = getRegionIJ(GEAR.X_coor,GEAR.Y_coor,region.minE,region.minN);
[i2,j2] = getRegionIJ(GEAR.X_coor,GEAR.Y_coor,region.minE+region.dx*(region.dimE-1),...
    region.minN+region.dx*(region.dimN-1));
GEARthis = GEAR.val(:,i1:i2,j1:j2);


year = 1980:1999;
sim = squeeze(nanmean(reshape(CPMthis,size(CPMthis,1),size(CPMthis,2),1,[]),4))';
plot(year,sim,'k-','linewidth',1);
hold on;
obs = squeeze(nanmean(reshape(GEARthis,size(GEARthis,1),1,[]),3));
plot(year,obs,'r-');
title(region.Name);

% % remove LOCs far from RGs
% [badI] = find(GEAR.dist(i1:i2,j1:j2)>2000);
% GEARthis(:,badI) = NaN;

GEARthis = imresize3(GEARthis,'Scale',[1,5,5],'method','nearest');
GEARthis = imresize3(GEARthis,'Scale',[1,1/11,1/11],'method','box');

STATS = table;
STATS.obs_cv = {squeeze(nanstd(GEARthis,0,1)./nanmean(GEARthis,1))};
STATS.sim_cv = {squeeze(nanstd(CPMthis,0,2)./nanmean(CPMthis,2))};

end


function plot_ENESEMBLE_MONTHLY_MEAN();
% ensembles monthly mean.

for mon = 1:12
    
    [CPM,E,N] = getMonthCPMEns(region,mon);% [mon,ensNo]
    
    
    setFigureProperty('Meeting');
    ha = tight_subplot(3,4,[.05 .05],[.05 .05],[.15 .05]);
    
    for EnsNo = 1:numel(ENSEMBLENO)
        
        axes(ha(EnsNo))
        in = getTrimTag('unit','km','product','cpm2.2');
        cpm = nanmean(cat(3,CPM{1,EnsNo}), 3)*24*365;
        cpm(~in) = NaN;
        
        
        plotAllAvera(mon,[],[],[],E,N,cpm,[])
        
        text(300,1250,sprintf('%s-Ens%s',getMonthName(mon),ENSEMBLENO{EnsNo}),...
            'fontsize',13,'horizontalalignment','center')
    end
    
    
    c = colorbar('location','Manual', 'position', [0.1 0.08 0.02 0.84]);
    
    c.Ruler.TickLabelFormat='%gmm';
    filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\ensemblePattern';
    mkdir(filePath);
    filename = [filePath,filesep,'cpm_ensemble1_12_Mon',sprintf('%02d',mon)];
    savePlot(filename,'XYWH',[150,100,700,600],'needreply','N');
    
    close all;
    
end

end



function plotAllAvera(mon,X_coor,Y_coor,gear1km,E,N,cpm,bias)


Z = cpm;
pcolor(E,N,Z);shading flat
cptcmap('precip_annualMeanUK', 'mapping','direct');%,'ncol',20);
xlim([0,700])
ylim([0,1200])
axis off

% axes(ha(2))
% pcolor(E,N,bias./Z*100);shading flat
% cptcmap('diff_darkBlue_darkRed', 'mapping','scaled','ncol',11);
% caxis([-100,100])
% xlim([0,700])
% ylim([0,1200])
% axis off
% text(200,1100,sprintf('Bias'))


UKMap = load('H:\CODE_MATLAB\SpatialTemporalDATA\shapeFileFolder\UKBorderGrid.mat');
hold on; plot(UKMap.borderE/1000,UKMap.borderN/1000,'k-');

end






