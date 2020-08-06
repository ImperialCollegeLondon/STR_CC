% ----------------------------------------------------------------------- %
% THIS SCRIPT IS TO VISUALIZE SEVERAL STATISTICS FOR CONVECTIVE STROM DURING
% SUMMERTIME
%
% This script is mainly to visualise wet spells & dry spells (pdf/cdf of each)
%
% @ Yuting Chen
% Imperial College London
% yuting.chen17@imperial.ac.uk
% update: 2020.08.06
% ----------------------------------------------------------------------- %

global regionName yi regCol
ENSEMBLENO = getEnsNos();
ha = tight_subplot(3,2,[.02 .1],[.1 .1],[.10 .05]);
itag = 1;
% Several Config
set(gcf,'units','points','position',[50,0,440,675]);
for regionName = {'CPM_NW','CPM_NE','CPM_S'}%{'EUK','WAL','SCO'}
    regionName = regionName{1};
    setFigureProperty('Paper');
        for y_name = {'wetSpell','drySpell'}
        y_name = y_name{1};
        axes(ha(itag))
        
        yi = 1:100; % yi = yi*100;
        for FP = {'2060-2080'}%{'2020-2040','2060-2080'}%{'2007-2018'}%
            FP = FP{1};
            
            [CPM,CPM_fut,RAD] = getAllforSpells(regionName,ENSEMBLENO,y_name,FP);
            
            [h,H_diff,P_diff] = plotOne(CPM,CPM_fut,RAD);
            text(35,0.08,[getLabel(y_name),'-',getPrintName(regionName)],'fontsize',10);
            if itag ~= 1
            legend off
            end
            % figure;boxplot(P_diff);median(P_diff)
            % title(regionName)
            %
            xlabel(getLabel(y_name));
            ylabel('PDF')
            itag = itag+1;
            % savePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\StormTracking';
            % fileName = sprintf('UKCP_plot_CStorm_JJA_%s_%s_dif%s_%s-%s',y_name,regionName,...
            %     FP,ENSEMBLENO{1},ENSEMBLENO{end});
            % savePlot([savePath,filesep,fileName],'XYWH',[50,0,220,200],'needreply','N','onlyPng',false);
            % pause(0.1)
            % close all
            
        end
    end
end

format_xylabel(ha,3,2)
savePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\StormTracking';
fileName = sprintf('UKCP_plot_CStorm_JJA_%s_dif%s_%s-%s',y_name,...
    FP,ENSEMBLENO{1},ENSEMBLENO{end});
savePlot([savePath,filesep,fileName],'XYWH',[50,0,500,675],'needreply','Y','onlyPng',false);

%% PLOT for other variables. % (not used)

for regionName = {'CPM_NW','CPM_NE','CPM_S'}%{'SCO','WAL','EUK'}
    regionName = regionName{1};
    setFigureProperty('Paper');
    % Several Config
    for y_name = {'rvol'}% {'rdur','rvol'};
        y_name = y_name{1};
        yi = [1:100]*1000;
        for FP = {'2060-2080'}%{'2020-2040'}%{'2007-2018'}
            FP = FP{1};
            
            [CPM,CPM_fut,RAD] = getAllforHist(regionName,ENSEMBLENO,y_name,FP);
            
            [h,H_diff,P_diff] = plotOne(CPM,CPM_fut,RAD);
            
            set(gcf,'units','points','position',[50,0,350,500]);
            % figure;boxplot(P_diff);median(P_diff)
            savePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\StormTracking';
            fileName = sprintf('UKCP_plot_CStorm_JJA_%s_%s_dif%s_%s-%s',y_name,regionName,...
                FP,ENSEMBLENO{1},ENSEMBLENO{end});
            title([y_name,'-',getPrintName(regionName)]);
            % savePlot([savePath,filesep,fileName],'XYWH',[50,0,350,500],'needreply','N','onlyPng',false);
            pause(0.1)
            
            hold on
            % close all
            
        end
    end
end



%% AUXILLARY FUNCTION
function [CPM,CPM_fut,RAD] = getAllforSpells(regionName,ENSEMBLENO,y_name,FP);
[STATS0,STATS] = deal([]);
for enInd = 1:length(ENSEMBLENO)
    ensNo = ENSEMBLENO{enInd};
    [STATS0] = [STATS0;get4Spells('1980-2000',ensNo)];
    [STATS] = [STATS;get4Spells(FP,ensNo)];
end

STATS0(STATS0.rpmax<5,:) = [];
STATS(STATS.rpmax<5,:) = [];

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

function [CPM,CPM_fut,RAD] = getAllforHist(regionName,ENSEMBLENO,y_name,FP);
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

function label = getLabel(x_name)
switch(x_name)
    case 'drySpell'
        label = 'Dry Spells (hour)';
    case 'wetSpell'
        label = 'Wet Spells (hour)';
    otherwise
        error('check Function getLabel(x_name)');
end

end


function [A] = get4Spells(Period,ensNo)
global regionName

STATS = [];
MON = 6:8;
Config = getConfig(upper(regionName),MON,Period,ensNo);
if strcmp(Period,'2007-2018')
    Dat = load([Config.saveIt.path,filesep,sprintf('CS_%s_%s_FIELD_%02d-%02d_%s.mat',...
        regionName,Period,Config.Month(1),Config.Month(end),'RAD')],...
        'TE','Config');
else
    Dat = load([Config.saveIt.path,filesep,sprintf('CS_%s_%s_FIELD_%02d-%02d_%s.mat',...
        regionName,Period,Config.Month(1),Config.Month(end),ensNo)],...
        'TE','Config');
end

sT = (cellfun(@(x)x(1),Dat.TE,'UniformOutput', false));
sT = cat(1,[],sT{:});
eT = (cellfun(@(x)x(end),Dat.TE,'UniformOutput', false));
eT = cat(1,[],eT{:});
wetSpell = hours(eT-sT)+1;
drySpell = [NaN;hours(sT(2:end)-eT(1:end-1))];
drySpell(drySpell>hours(days(6*30))) = NaN;
wetSpell(wetSpell>hours(days(6*30))) = NaN;

A = table(wetSpell,drySpell);

A.ensNo = repmat(string(ensNo),size(A.wetSpell));

[B] = get4Plot_all(Period,ensNo);
A.rpmax = B.rpmax;
A.rdur = B.rdur;

    function [A] = get4Plot_all(Period,ensNo)
        
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
        
        % A = STATS(STATS.rpmax>=5 & ~isnan(STATS.rspeed),:);
        % A = STATS(STATS.rpmax>=5,:);
        
        A = table();
        STATS.tag = STATS.evi;
        STATS.mon = [];
        D = grpstats(STATS,{'tag'},{'mean','max','median'});
        A.rpmax = D.max_rpmax;
        A.rdur = D.median_rdur;
        
    end

end


function A = get4Plot(Period,ensNo)
global regionName

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

A = STATS(STATS.rsize>=2.2*2.2 & STATS.rpmax>5 & ~isnan(STATS.rspeed),:);

% A = table();
% STATS.tag = STATS.evi;
% STATS.mon = [];
% D = grpstats(STATS,{'tag'},{'mean','max','median'});
% A.rvol = D.max_rvol; % mean_rvol.*D.median_rdur/24;%
% A.rsize = D.max_rsize;
% A.rpmax = D.max_rpmax;
% A.rspeed = D.median_rspeed;
% A.rdur = D.median_rdur;
% A(A.rpmax<5,:) = [];

A.ensNo = repmat(string(ensNo),size(A.rpmax));

end


function [ax,H_diff,P_diff] = plotOne(CPM,CPM_fut,RAD)

global yi regionName regCol
bw = 3;
sm = @(x)x;%reshape(smooth(x,bw),size(x));
kernelName = 'normal';%'box';%
obsCol = [66 146 199]/255;
futCol = [240  60  43]/255;
f = [];
for ensNo = 1:12
    A = CPM{ensNo};
    [N,EDGES] = histcounts(A,[yi,inf],'Normalization','pdf');N = conv(N,ones(1,5)/5,'same');
    f(ensNo,:) = N;
end
f = f+1e-8;
ax = gca;
plot(yi,sm(nanmedian(f,1)),'-','color',obsCol);
hold on;
ylim([0,0.25])%%%%%%%%%%%%%%

hsimran = fill([yi,flip(yi)],[sm(nanmin(f,[],1)),sm(flip(nanmax(f,[],1)))],...
    obsCol,'LineStyle','none');%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f_fut = [];
for ensNo = 1:12
    A = CPM_fut{ensNo};
    [N,EDGES] = histcounts(A,[yi,inf],'Normalization','pdf');N = conv(N,ones(1,5)/5,'same');
    f_fut(ensNo,:) = N;
end
f_fut = f_fut+1e-8;
plot(yi,sm(nanmedian(f_fut,1)),'-','color',futCol);
hfutran = fill([yi,flip(yi)],[sm(nanmin(f_fut,[],1)),sm(flip(nanmax(f_fut,[],1)))],futCol,...
    'LineStyle','none');

alpha(0.3)

[H_diff,P_diff] = deal([]);
for ensNo = 1:12
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
set(gca,'YScale','log')
ax = gca;
ax.YLim = [4e-4,0.1];
ax.YTick = [10^-3,10^-2,10^-1];
ax.XTick = [1,24,48,72,96];
hold off;
grid on;
axis('square')
legend([hsimran,hfutran],strcat('past-','(1980-2000)'),strcat('future-','(2060-2080)'))
legend boxoff
end

function pn = getPrintName(regionName)
switch(regionName)
    case 'WAL'
        pn = 'SWUK';
    case 'SCO'
        pn = 'NUK';
    case 'EUK'
        pn = 'SEUK';
    case 'CPM_NW'
        pn = 'NW-UK';
    case 'CPM_NE'
        pn = 'NE-UK';
    case 'CPM_S'
        pn = 'S-UK';
end
end