% ----------------------------------------------------------------------- %
% THIS SCRIPT IS TO PLOT QUANTILE DIFFERENCES OF SEVERAL STORM PROPERTIES
% (HOURLY SCALE)
%
% @ Yuting Chen
% Imperial College London
% yuting.chen17@imperial.ac.uk
% UPDATE: 2020.08.06
% ----------------------------------------------------------------------- %

close all
global regionName yi in
ENSEMBLENO = getEnsNos();
itag = 1;
for y_name = {'rpmax','rsize','rami','rspeed'}%{'rrmi','rsizeall',}% 
    subplot(2,2,itag)
    y_name = y_name{1};
    legh = handle(0);
    for regionName = {'CPM_NW','CPM_NE','CPM_S'}%{'SCO','WAL','EUK'}
        regionName = regionName{1};
        setFigureProperty('Paper');
        
        [E1,N1] = getEN(getfield(REGIONS_info(),regionName));
        UKMap = getUKMap();
        in = inpolygon(E1,N1,UKMap.borderE/1000,UKMap.borderN/1000);
        
        % Several Config
        
        for FP = {'2060-2080'}%{'2020-2040','2060-2080'}%{'2007-2018'}%
            FP = FP{1};
            
            [CPM,CPM_fut,RAD] = getDataforHist(regionName,ENSEMBLENO,y_name,FP);
            
            [h,H_diff,P_diff,legh(end+1)] = plotOne(CPM,CPM_fut,RAD);
            
            pause(0.1)
            grid off
            box off
            hold on
            ax = gca;
            ax.YLim = [-0.25,0.5]*100;
            ax.YTick = -20:10:50;
        end
    end
    legh(1) = [];
    
    set(gcf,'units','centimeters','position',[0 0 8 8]);
    plot(ax.XLim,[0,0],'-','color',[0 0 0 0.5],'linewidth',0.5);
    legend(legh,{'NW-UK','NE-UK','S-UK'},'Location','NorthWest')
    legend boxoff
    axis('square')
    savePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\StormTracking';
    fileName = sprintf('UKCP_plot_CStorm_JJA_%s_%s_dif%s_%s-%s','cpmRegions',y_name,...
        FP,ENSEMBLENO{1},ENSEMBLENO{end});
    xlabel(getLabel(y_name))
    ylabel('quantile difference[%]')
    ytickformat(ax, 'percentage');
    
    itag = itag+1;
end

fileName = sprintf('UKCP_plot_CStorm_JJA_%s_dif%s_%s-%s','cpmRegions',...
        FP,ENSEMBLENO{1},ENSEMBLENO{end});
savePlot([savePath,filesep,fileName],'units','centimeters','XYWH',[0 0 8*2.5 8*2.5],...
    'needreply','Y','onlyPng',false);


%% AUXILLARY FUNCTION
function [CPM,CPM_fut,RAD] = getDataforHist(regionName,ENSEMBLENO,y_name,FP);
[STATS0,STATS] = deal([]);
for enInd = 1:length(ENSEMBLENO)
    ensNo = ENSEMBLENO{enInd};
    STATS0 = [STATS0;get4Plot('1980-2000',ensNo)];
    STATS = [STATS;get4Plot(FP,ensNo)];
end

[CPM,CPM_fut] = deal([]);
for enInd = 1:length(ENSEMBLENO)
    ensNo = ENSEMBLENO{enInd};
    
    thisSTATS = STATS0(strcmp(STATS0.ensNo,ensNo),:);
    CPM{enInd} = getfield(thisSTATS,y_name);
    
    thisSTATS = STATS(strcmp(STATS.ensNo,ensNo),:);
    CPM_fut{enInd} = getfield(thisSTATS,y_name);
    
    RAD = [];
end
end

function A = get4Plot(Period,ensNo)
global regionName in

STATS = [];
MON = 6:8;
Config = getConfig(upper(regionName),MON,Period,ensNo);
if strcmp(Period,'2007-2018')
    Dat = load([Config.saveIt.path,filesep,sprintf('CS_%s_%s_STATS_%02d-%02d_%s.mat',...
        regionName,Period,Config.Month(1),Config.Month(end),'RAD')],...
        'STATS','Config');
else
    Dat = load([Config.saveIt.path,filesep,sprintf('CS_%s_%s_STATS_%02d-%02d_%s.mat',...
        regionName,Period,Config.Month(1),Config.Month(end),ensNo)],...
        'STATS','Config');
end
STATS = [STATS;Dat.STATS];

A = STATS(STATS.rpmax>=10,:);% STATS(STATS.rpmax>=10 & STATS.rsize>5,:);
A.rami = A.rvol /getSumArea(regionName,in) /1000 *3600; % mm/h

A.ensNo = repmat(string(ensNo),size(A.rpmax));

end

function [ax,H_diff,P_diff,hleg] = plotOne(CPM,CPM_fut,RAD)

global yi regionName
regCol = getRegionColo(regionName);

sm = @(x)x;%reshape(smooth(x,bw),size(x));
kernelName = 'normal';
hold on;
f = [];
yi = [];
A = [];
for ensNo = 1:length(CPM)
    A = cat(1,A,CPM{ensNo});
end
quantiles = [0.01:0.01:0.99];% [0.01:0.01:0.99];
yi = quantile(A,quantiles);
for ensNo = 1:length(CPM)
    A = CPM{ensNo};
    f(ensNo,:) = quantile(A,quantiles);
end
ax = gca;

f_fut = [];
for ensNo = 1:length(CPM)
    A = CPM_fut{ensNo};
    f_fut(ensNo,:) = quantile(A,quantiles);
end

relDif = 100*(f_fut-f)./f;
hleg = plot(yi,sm(nanmedian(relDif,1)),'-','marker',getMarker(regionName),'markerfacecolor',[0.9,0.9,0.9],...
    'color',[regCol,0.8],'linewidth',0.5,'markersize',3);
for ii = 1:length(yi)
    hfutran = plot([yi(ii),yi(ii)],[sm(prctile(relDif(:,ii),25)),sm(prctile(relDif(:,ii),75))],'-',...
        'Color',[regCol,0.5],'linewidth',0.1);
end

alpha(0.1)

[H_diff,P_diff] = deal([]);
for ensNo = 1:length(CPM)
    [H_diff(ensNo),P_diff(ensNo)] = kstest2(CPM{ensNo},CPM_fut{ensNo},...
        'tail','unequal');
    [P_diff(ensNo),H_diff(ensNo)] = ranksum(CPM{ensNo},CPM_fut{ensNo},...
        'tail','right');
end


if ~isempty(RAD)
    A = RAD;
    [N,EDGES] = histcounts(A,[yi,inf],'Normalization','pdf');
    f(ensNo,:) = N;
    plot(yi,sm(fobs),'k-');
end
box on;
ax = gca;
% ax.YLim = [4e-4,0.1];
hold off;
grid on;

end
function mStr = getMarker(regionName)
switch(regionName)
    case 'SCO'
        mStr = 's';
    case 'WAL'
        mStr = '^';
    case 'EUK'
        mStr = '>';
    case 'CPM_NW'
        mStr = 's';
    case 'CPM_NE'
        mStr = '^';
    case 'CPM_S'
        mStr = '>';
end
end
function regCol = getRegionColo(regionName)
switch(regionName)
    case 'SCO'
        regCol = [66 146 199]/255;
    case 'WAL'
        regCol = [240  60  43]/255;
    case 'EUK'
        regCol = [0.3 0.3 0.3];
    case 'CPM_NW'
        regCol = [66 146 199]/255;
    case 'CPM_NE'
        regCol = [240  60  43]/255;
    case 'CPM_S'
        regCol = [0.3 0.3 0.3];
end
end
function label = getLabel(x_name)
switch(x_name)
    case 'rdur'
        label = 'Duration(h)';
    case 'rspeed'
        label = 'Speed(km h^{-1})';
    case 'rvol'
        label = 'Ptotal(m^3 s^{-1})';
    case 'rrmi'
        label = 'rmi (mm)';
    case 'rpmax'
        label = 'Pmax(mm h^{-1})';
    case 'rsize'
        label = 'Size(km^2)';
    case 'rsizeall'
        label = 'SizeAll(km^2)';
    case 'rami'
        label = 'Areal Mean Intensity(mm/h)';
    otherwise
        error('check Function getLabel(x_name)');
end

end




