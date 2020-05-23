% This file is to map annual mean precipitation from UKCP CPM datasets
%
% Yuting Chen
%
% yuting.chen17@imperial.ac.uk
% Imperial College London



clear;clc

REGIONS = REGIONS_info();
region = REGIONS.UK;%Scotland;% Westuk;% wales

mons = 1:12;%6:8;

% Result from CPM
[CPM,E,N] = getMonthCPM(region);% [mon,ensNo]


in = getTrimTag('unit','km','product','cpm2.2');
cpm = 0;
for ensNo = 1:12
    cpm0 = nansum(cat(3,CPM{mons,ensNo}).*reshape(eomday(1999,mons),...
        [1,1,length(mons)]), 3)*24;
    cpm0(~in) = NaN;
    cpm = cpm0+cpm;
end
cpm = cpm/12;


sf = 365/sum(eomday(1999,mons));

%%
figure;
hold on;

UKMap = getUKMap();

set(gca,'Linewidth',2);
Z = cpm*sf;
pcolor(E,N,Z);shading flat
cptcmap('precip_annualMeanUK', 'mapping','direct');%,'ncol',20);

axis off
c = colorbar('location','Manual', 'position', [0.65 0.45 0.02 0.3],'fontsize',12);

c.Ruler.TickLabelFormat='%gmm';
c.FontSize = 12;
plot(UKMap.borderE/1000,UKMap.borderN/1000,'-','linewidth',1,'color',[0.5,0.5,0.5]);
ylim([0,1200]);caxis([0,3100]);xlim([0,800]);

REGIONS = REGIONS_info();
% plot(borderE/1000,borderN/1000,'k-','linewidth',0.5)
plotRec = @(region)rectangle('Position',[region.minE,region.minN,...
    region.dimE*region.dx,region.dimN*region.dx],'linewidth',1,'Linestyle','--');
plotRec(REGIONS.SCO);
plotRec(REGIONS.WAL);
plotRec(REGIONS.EUK);

text(270,770,'NUK','fontsize',12,'horizontalalignment','center','fontweight','bold','background',[0.9 0.9 0.9 0.6]);
text(204,220,'SWUK','fontsize',12,'horizontalalignment','center','fontweight','bold','background',[0.9 0.9 0.9 0.6]);
text(615,165,'SEUK','fontsize',12,'horizontalalignment','center','fontweight','bold','background',[0.9 0.9 0.9 0.6]);
text(200,-50,{'Background:CPM(1980-2000)';'AnnualMean'},'fontsize',10,'fontangle','italic')

hold off
axis equal
xlim([0,800]);

savePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\';
fileName = sprintf('StudyRegion');
savePlot([savePath,filesep,fileName],'units','centimeters','XYWH',[5,0,8,11],'needreply','Y','onlyPng',false);


