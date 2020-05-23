dataSP = 'H:\DATA_CLIMATE\UKCP18\';
REGIONS = REGIONS_info();


Name = {'SCO','WAL','EUK'};

figure;
rowno = 1;
% # separate plot # 
% colno = numel(Name);
% ha = tight_subplot(rowno,colno,[.08 .02],[.15 .05],[.15 .05]);
ha = figure; colno = 1;
setFigureProperty('Paper')
set(0,'defaultAxesFontSize',14,'defaultAxesFontName','Arial',...
    'defaultAxesTitleFontWeight','Bold',...
    'defaultTextFontSize',10)

global seasons yi
seasons = 2;% 1:4;% 1:4;% 1:4
hsim = handle(1);
for thisstation = 1:length(Name)
    
    % [CPM_pas,CPM_fut] = loadData(dataSp,Name{thisstation});
    [CPM_pas] = loadData(Name{thisstation},'1980-2000');
    [CPM_fut] = loadData(Name{thisstation},'2060-2080');
    
    minIMF = 5;
    maxIMF = 10000;
    % 'rvol','rsize','rdur'
%     % rpmax
%     yi = 5:0.5:110;%5:1:100;%
%     % axes(ha(thisstation+length(Name)*0))
%     [ax,cpmPasNum,cpmFutNum,hsim(thisstation)] = plotOne(CPM_pas,CPM_fut,'rpmax',minIMF,maxIMF);
    
    rpmaxVec = [5,40,50,60,70,80,90];
    [XVal,xlabels] = deal([]);
    for bi = 1:length(rpmaxVec)-1
        minP = rpmaxVec(bi);
        maxP = rpmaxVec(bi+1);
        if bi+1 == length(rpmaxVec)
            maxP = Inf;
        end
        [ax,cpmPasNum,cpmFutNum,~] = plotOne(CPM_pas,CPM_fut,'rpmax',minP,maxP);
        XVal(:,bi) = nanmedian((cpmFutNum-cpmPasNum)./cpmPasNum);
        if bi == 1
            xlabels{bi} = sprintf('<%d',maxP);
        elseif bi == length(rpmaxVec)-1
            xlabels{bi} = sprintf('>%d',minP);
        else
            xlabels{bi} = sprintf('%d-%d',minP,maxP);
        end
    end
    sm = @(x)x;%reshape(smooth(x,bw),size(x));
    
    hsim(thisstation) = plot(rpmaxVec(2:end),...
        100*sm(nanmedian(XVal,1)),'-');% ylim([0,0.25])
    hold on;
    plot(ax.XLim,[0,0],'k--','Linewidth',0.5);
    hold on
    ax.XTick= [rpmaxVec(2:end)];
    ax.XTickLabel = xlabels;
    ylim([-50,400]) % in percentage
    xlim([35,95])
    
    % ax.YLim=[1e-5,0.07];%[0,1];%
    formatIt('PMax [mm/h]',Name{thisstation})
    % ax.XTick= [5,25,50,75,100];
    % title(Name{thisstation})
    
    fprintf('---Occurrence of rainfall between [%3d,%3d]mm/h---\n',minIMF,maxIMF);
    fprintf('Region %s Past:%d/yr \t Future:%d/yr \t Change:%.3f \n',...
        Name{thisstation},round(nanmedian(cpmPasNum/20)),...
        round(nanmedian(cpmFutNum/20)),...
        nanmedian((cpmFutNum-cpmPasNum)./cpmPasNum));
    
    
end
xtickangle(90)
ytickformat('percentage');
legend([hsim],Name,'Location','Northwest')
% legend([hsimran,hsimfutran,hobs],{'CPM(1980-2000)','CPM(2060-2080)','Radar(2007-2018)'})

% legend([hsimran,hobs],{'CPM(1980-2000)','Radar(2007-2018)'})

set(ha(((1:rowno)-1)'*colno+(2:colno)),'YTickLabel',[],'YLabel',[])

filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\';
filename = [filePath,filesep,'RChangePDFs_event_',getSeasonName(seasons(1)),'to',getSeasonName(seasons(end))];
savePlot(filename,'units','centimeters','XYWH',[5,0,9,9],'needreply','Y');


function [CPM] = loadData(regionName,Period)
CPM = [];
MON = 6:8;
for ensNo = getEnsNos()
    ensNo = ensNo{1};
    STATS = [];
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
    
    
    A = STATS(STATS.rpmax>=5,:);
%     A = STATS(STATS.rpmax>=5 & ~isnan(STATS.rspeed),:);
    
%     A = table();
%     STATS.tag = STATS.evi;
%     STATS.mon = [];
%     D = grpstats(STATS,{'tag'},{'mean','max','median'});
%     A.rvol = D.max_rvol; % mean_rvol.*D.median_rdur/24;%
%     A.rsize = D.max_rsize;
%     A.rpmax = D.max_rpmax;
%     A.rspeed = D.median_rspeed;
%     A.rdur = D.median_rdur;
%     A.ensNo = repmat(string(ensNo),size(A.rpmax));
%     STATS = A;
    
    CPM = cat(2,CPM,{STATS});  
end
end

function formatIt(xlab,name0)
ylabel('Relative Change of PDF');
set(gca,'linewidth',2)
xlabel(xlab)
ax = gca;
% text(ax.XLim(2)*0.9,ax.YLim(2),name0,'HorizontalAlignment',...
%     'Right','VerticalAlignment','Top');
end

function titleIt(name0)
if strcmp(name0,'London')
    name = 'SouthUK';
end
title(sprintf('%s',name0))% (110x110Km^2)
end



function [ax,cpmPasNum,cpmFutNum,hsim] = plotOne(CPM_pas,CPM_fut,var,minIMF,maxIMF)
% Input
%      var: 'rvol','rsize','rdur','rpmax'
%

global seasons yi
hsim = [];
bw = 10;%%(yi(2)-yi(1))*10;
if yi(2)-yi(1) == 1
    bw = 1;
end
sm = @(x)x;%reshape(smooth(x,bw),size(x));
kernelName = 'normal';%'box';%

% hold on;

[f_pas,~,cpmPasNum] = plotEnsPdf(CPM_pas);
[f_fut,~,cpmFutNum] = plotEnsPdf(CPM_fut);
ax = gca;

% plot(yi,sm(nanmedian(f_pas,1)),'b-');% ylim([0,0.25])
% hsimpas = fill([yi,flip(yi)],[sm(nanmin(f_pas,[],1)),sm(flip(nanmax(f_pas,[],1)))],'b',...
%     'LineStyle','none');
% plot(yi,sm(nanmedian(f_pas,1)),'r-');% ylim([0,0.25])
% hsimfut = fill([yi,flip(yi)],[sm(nanmin(f_pas,[],1)),sm(flip(nanmax(f_pas,[],1)))],'r',...
%     'LineStyle','none');
%  
% alpha(0.3)
% box on;
% hold off;

    function [f,xi,cpmNum] =  plotEnsPdf(CPM_ENS)
        [f,xi] = deal([]);
        cpmNum = [];
        for ensNo = 1:12
            A = [];
            A0 = getfield(CPM_ENS{ensNo},var);
            rmind = CPM_ENS{ensNo}.rpmax<minIMF | CPM_ENS{ensNo}.rpmax>maxIMF | A0 == 0;
            A0(rmind) = [];
            A = cat(1,A,A0);
            cpmNum(ensNo) = numel(A);
            if yi(2)-yi(1) == 1
                [f(ensNo,:),xi(ensNo,:)] = ksdensity(A,yi,'function','pdf','kernel',kernelName,...
                    'support', 'positive', 'BoundaryCorrection',...
                    'Reflection','Bandwidth',2);
            else
                [f(ensNo,:),xi(ensNo,:),bw0]= ksdensity(A,yi,'function','pdf','kernel',kernelName,...
                    'support', 'positive', 'BoundaryCorrection',...
                    'Reflection');
                [f(ensNo,:),xi(ensNo,:),~]= ksdensity(A,yi,'function','pdf','kernel',kernelName,...
                    'support', 'positive', 'BoundaryCorrection',...
                    'Reflection','Bandwidth',10);%bw0/2);
            end
        end
    end



end


