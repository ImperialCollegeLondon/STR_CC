clear;clc
close all

durVec = [1 3 6 12 24];

for regionName = {'SCO','EUK','WAL'}%{'WAL'}% 
    
    regionName = regionName{1};
    
    dataSourcePath = 'D:\UKCP18\fixedAreaARF_oldRegions';%'D:\UKCP18\stormcentredARF';%
    close all
    figure;
    setFigureProperty('Subplot2')
    
   %  for ensNo = getEnsNos();
   % [Dat_past,Dat_future,arf] = getData(dataSourcePath,regionName,ensNo);
    [Dat_past,Dat_radar,arf] = getData(dataSourcePath,regionName,getEnsNos());

    subplot 121
    hold on;
    [hpas,hobs] = deal(handle(0));
    for duri = 1:length(durVec)
        dur = durVec(duri);
        colobs = getColorMat(duri,'rad');
        colpas = getColorMat(duri,'pas');
        % colfut = getColorMat(duri,'fut');
        val = Dat_past.arf(Dat_past.dur==dur,:);
        hpas(duri) = plot([1,arf.area]*2.2*2.2,[1,nanmean(val,1)],...
            '-','color',[colpas,0.5],'markerfacecolor',[colpas],'markersize',5);
        val = Dat_radar.arf(Dat_radar.dur==dur,:);
        hobs(duri) = plot([1,arf.area]*2.2*2.2,[1,nanmean(val,1)],...
            '--.','color',[colobs,0.5],'markerfacecolor',[colobs],'markersize',10);
    end
    formatSubP1(hpas(1),hobs(1),durVec([1]))
    hold off
   
    subplot 122
    hold on;
    [hpas,hobs] = deal(handle(0));
    for duri = 1:length(durVec)
        dur = durVec(duri);
        col = getColorMat(duri,'rad');
        hpas(duri) = plot([arf.area]*2.2*2.2,...
            nanmean(Dat_radar.arf(Dat_radar.dur==dur,:),1)./nanmean(Dat_past.arf(Dat_past.dur==dur,:),1),...
            'o-','color',[col,0.5],'markerfacecolor',[col]);
    end
    formatSubP2(hpas(1),hobs(1),durVec([1]))
    hold off
    
    figPath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\ARF';
    fileName = sprintf('%s%sARF_NERC_JJA_withObs_%s',figPath,filesep,regionName);
    savePlot(fileName,'XYWH',[150,0,550,250],'needreply','N');
    
end
%% AUXILLARY FUNCTION
function colorMat = getColorMat(durVec,type)
switch(type)
    case 'pas'
        colorMat = [66 146 199]/255;
    case 'fut'
        colorMat = [240  60  43]/255;
    case 'rad'
        colorMat = [0.5 0.5 0.5];
end
end

function [A,B,arf] = getData(dataSourcePath,regionName,ENSEMBLENO)
region = getfield(REGIONS_info(),regionName);
[A] = deal([]);

% cpm
for ensNo = ENSEMBLENO
    ensNo = ensNo{1};
    try
        fileName = sprintf('%s%sARFs_%s_ensNo%s_%s.mat',dataSourcePath,filesep,...
            region.Name,ensNo,'1980-2000');
        load(fileName,'T');
        T.allT_Max = [];
        T(~any(~isnan(T.arf'))',:) = [];
        A = [A;T];
    catch me
        break
    end
end

% radar
fileName = sprintf('%s%sARFs_%s_ensNo%s_%s.mat',dataSourcePath,filesep,...
    region.Name,'RAD','2007-2018');
load(fileName,'T');
T.allT_Max = [];
T(~any(~isnan(T.arf'))',:) = [];
B = T;

arf = struct;
arf.area = [10,30,80,150,250,380,530,700,900,1250,2000,4000,8000]./(2.2).^2;
% arf.area = [10,30,80,150,250,380,530,700,900]./(2.2).^2;

end
function formatSubP2(hpas,hobs,durVec)
ylim([0.85,1.15])
grid minor
box on;
ax = gca;
ax.LineWidth = 1;
ax.XScale = 'log';
ax.XTick = [10,100,1000,8000];
ax.XLim = [10,5000];
plot(ax.XLim,[1,1],'k--','color',[0.1,0.1,0.1,0.5]);
legend(flip([hpas]),flip(strcat(string(durVec'),'h')),...
    'Location','SouthWest')
legend boxoff
xlabel('Area (km^2)')
ylabel('radar/cpm2.2');
end

function formatSubP1(hpas,hobs,durVec)
ylim([0.5,1])
grid minor
box on;
ax = gca;
legend(flip([hpas,hobs]),flip({strcat(string(durVec'),'h-cpm2.2'),...
    strcat(string(durVec'),'h-radar')}),...
    'Location','SouthWest')
legend boxoff
xlabel('Area (km^2)')
ylabel('Areal Reduction Factor');
ax.XScale = 'log';
ax.XTick = [10,100,1000,8000];
ax.XLim = [10,5000];
ax.LineWidth = 1;
end

