% ---------------------------------------------------------------------- %
% This file achieves several figures which are used in D3.3 Report
% MAINLY INCLUDE THOSE FIGURES IN CHAPTER OF MODEL CONFIGURATION
% @ Yuting Chen
% yuting.chen17@imperial.ac.uk
% Imperial College London
% ----------------------------------------------------------------------- %


clear;clc

datSavePath = 'K:\DATA_FCS\CrossVal';
figSavePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_FCS';
% % 'C2P5-1',
for expNo = {'C1P5-Rad1','C1P5In-Rad1','C1P5In-1','C1P5-1',...
        'C1P5-2','C1U5-1','C1U10-1',...
        'C1U5-Rad1','C1U10-Rad1','C2P5-Rad1'}
    saveSummaryStats(expNo,datSavePath)
end

%% OUTPUT TABLE FOR COMPARISON WITH OTHER PILOTS

EXPS = {'C1P5RadComb'};%{'C1P5-Rad1'};% {'C1P5RadComb'};
[ETab,TP,FN,FP,TN,evi] = plotEvalMat_Mean(EXPS);

TP = round(TP);
FN = round(FN);
FP = round(FP);
TN = round(TN);
prc25 = @(x)prctile(x,25);
        prc75 = @(x)prctile(x,75);
getM = @(x)nanmedian(x);%nanmedian(x);
        
statarray = grpstats(table(TP,FN,FP,TN,evi),'evi',{'median','mode','std',@(x)prc25(x),@(x)prc75(x), @(x)std(x)./nanmean(x)});
%%
ETab_f = [getM(statarray.Fun4_TP),getM(statarray.median_TP),getM(statarray.Fun5_TP),getM(statarray.Fun6_TP);...
    getM(statarray.Fun4_FN),getM(statarray.median_FN),getM(statarray.Fun5_FN),getM(statarray.Fun6_FN);...
    getM(statarray.Fun4_FP),getM(statarray.median_FP),getM(statarray.Fun5_FP),getM(statarray.Fun6_FP);...
    getM(statarray.Fun4_TN),getM(statarray.median_TN),getM(statarray.Fun5_TN),getM(statarray.Fun6_TN)];

%%
mi = 2;
[TP,FN,FP,TN] = SolveIt(ETab(1,mi),ETab(2,mi),...
    ETab(3,mi),ETab(4,mi),ETab(5,mi));


%% LOOK AT TABLE SUMMARY

f = figure;
f.GraphicsSmoothing = 'on';
EXPS = {'C1P5RadComb'};%{'C1P5-Rad1'};% {'C1P5RadComb'};
ETab = plotEvalMat_50Prc(EXPS);


%% SPECIFIC EXPERIMENT + ALL CHECKING

EXPS = {'C1P5In-Rad1','C1P5In-1','C1P5-Rad1','C1P5-1'};

setFigureProperty('Subplot2');
plotEvalMat_50Prc(EXPS)

filename = [figSavePath,filesep,'ConfigurationOptimization_RainfallProduct_Radar_KED_Ind'];
savePlot(filename,'XYWH',[150,0,350,300],'needreply','N');


%% ALL EXPERIMENT + CHECKING ONE
% 'C1P5RadComb'
EXPS = {'C1P5-Rad1','C1P5In-Rad1','C1P5-1','C1P5In-1','C1U5-1','C1U5-Rad1','C1U10-1','C1U10-Rad1','C1P5-2'};
% MatGroup = {'Hit Rate','False Alarm','Critical Success Index','Error Bias'};
MatGroup = {'True positive rate','True negative rate',...
    'Positive predictive rate','Negative predictive rate',...
    'Accuracy'};
ALLMA = {'tpr','tnr','ppr','npr','acc'};

for i = 1:5
    evalM = ALLMA{i};
    
    figure;
    % [~,I] = ismember(evalM,{'hit','fa','cci','eb'});
    [~,I] = ismember(evalM,{'tpr','tnr','ppr','npr','acc'});
    
    set(gcf,'position',[150,0,200,300]);

    valthis = plotEvalOne(EXPS,evalM,MatGroup{I});

    filename = [figSavePath,filesep,'ConfigurationOptimization_AllExps_',evalM];
    savePlot(filename,'Units','centimeters','XYWH',[0,0,4.5,6],'needreply','N','onlyPng',true);
    close all
    
end

%% INVESTIGATE: INFLUENCE OF $TIME TO STORM ORIGIN$ ON $PREDICTION PERFORMANCE$

figure;
rn=1; cn=3;
% ha = tight_subplot(rn,cn,[0.0 0.0],[.15 .05],[.15 .05]);
setFigureProperty('Subplot2');
datPath = ['K:\DATA_FCS\CrossVal\'];
unitb = 10;
htag = 1;h = handle(0);hname = [];

MASTR =  {'o','^','-'};%{'o','^','s','d','>'};%
[colorm, ~, ~, ~, ~] = cptcmap('cw1-002','mapping','scaled','ncol',length(MASTR)+1);
% for expNo = {'C1P5-Rad1','C1P5In-Rad1','C1U10-Rad1','C1U5-Rad1'} % Perf(event duration)
for expNo = {'C1P5-Rad1','C1P5In-Rad1','C1P5RadComb'} % Perf(time to storm origin)
    
    expNo = expNo{1};
    makerstr = MASTR{htag};
    load([datPath,sprintf('SummarySTATS_%s_IDW%g_uniqueEvent%g_includePrevious%g.mat',...
        expNo,0,0,0)],...
        'STATS','testConfig');
    
    x = STATS.eventTT/12;
    xgroup = round(x/6)*6;% look at event duration
    
    x = [];
    xgroup = 100/unitb *round(unitb*STATS.eventTi./STATS.eventTT)/100;
%     
    
    cm = colorm(htag,:);
    
    %     subplot(rn,cn,1);ylim([0,40]);plotg(STATS.rmse,xgroup,'RMSE',makerstr);
    %
    %     subplot(rn,cn,2);plotg(STATS.mape,xgroup,'MAPE',makerstr);
    
    subplot(rn,cn,1)
    ylim([0,1]);h(htag) = plotg(STATS.tpr,x,xgroup,'TPR',makerstr,cm);
    % ax = gca;ax.XTick = [];
    
    subplot(rn,cn,2)
    ylim([0,1]);plotg(STATS.ppr,x,xgroup,'PPR',makerstr,cm);
    % ax = gca;ax.XTick = [];
    
    subplot(rn,cn,3)
    ylim([0,1]);plotg(STATS.acc,x,xgroup,'ACC',makerstr,cm);
    
%     subplot(rn,cn,4)
%     ylim([0,1]);plotg(STATS.npr,x,xgroup,'NPR',makerstr,cm);
%     % set(gca,'YScale','log')
%     subplot(rn,cn,5)
%     ylim([0,1]);plotg(STATS.tnr,x,xgroup,'TNR',makerstr,cm);
%     % set(gca,'YScale','log')

    hname{htag} = expNo;
    htag = htag+1;
end
legend(h,hname,'Location','SouthWest','Fontsize',8)

% format_xylabel(ha,rn,cn)
filename = [figSavePath,filesep,'offlineVali_combineC1P5_-Rad1'];
savePlot(filename,'Units','centimeters','XYWH',[0,0,4.5*3+4,5.5],'needreply','Y');
% close all
%

% filename = [figSavePath,filesep,'offlineVali_withD2StormOrigin_UMPCA'];
% savePlot(filename,'Units','centimeters','XYWH',[0,0,4.5*3+4,10],'needreply','Y');


%%
excluP = 0;
figure;
setFigureProperty('Subplot2');

ha = tight_subplot(4,1,[.08 .08],[.05 .05],[.05 .05]);

axes(ha(1))

formatHist2(STATS.hit(STATS.eventTi>excluP*STATS.eventTT),...
    STATS.eventNo(STATS.eventTi>excluP*STATS.eventTT),...
    'HIT RATE',[0,0.5,1]);ylim([0,1])

axes(ha(2))

formatHist2(STATS.fa(STATS.eventTi>excluP*STATS.eventTT),...
    STATS.eventNo(STATS.eventTi>excluP*STATS.eventTT),...
    'FALSE ALARM',[0,0.5,1]);ylim([0,1])

axes(ha(3))

formatHist2(STATS.cci(STATS.eventTi>excluP*STATS.eventTT),...
    STATS.eventNo(STATS.eventTi>excluP*STATS.eventTT),...
    'CRITICAL SUCCESS INDEX',[0,0.5,1]);ylim([0,1])

axes(ha(4))

formatHist2(STATS.eb(STATS.eventTi>excluP*STATS.eventTT),...
    STATS.eventNo(STATS.eventTi>excluP*STATS.eventTT),...
    'ERROR BIAS',[0,1,1e8]);ylim([0,2])

set(gcf,'units','points','position',[0 0 800 600]);
filename = [figSavePath,filesep,'offlineVali_HistStormOrigin_UMPCA'];
savePlot(filename,'XYWH',[0 0 800 600],'needreply','Y');


%% AUXILLARY FUNCTION
function h = plotg(Y,X,groupx,ylabelText,makerstr,cm)

hold on
y = grpstats(Y,groupx,'median');
h = plot(unique(groupx(:)),y,makerstr,'color',cm*0.8,'linewidth',1,'markerfacecolor',cm,...
    'markersize',4);
hold on;
if ~isempty(X)
    % fit poly1 to datasets.
    cfun = fit(X,Y,'poly1','Normalize','on','Robust','Bisquare');
    y = feval(cfun,groupx(:));
    plot(groupx(:),y,'--','color',cm,'linewidth',1);
end
% plot(unique(groupx(:)),grpstats(Y,groupx,@(x)prctile(x,75)),'--','linewidth',0.5);
% plot(unique(groupx(:)),grpstats(Y,groupx,@(x)prctile(x,25)),'--','linewidth',0.5);
if 0
    boxplot(Y,groupx,'whisker',100,'colors',copper(20))
    af = gcf;
    set(findobj(af,'LineStyle','--'),'LineStyle','-');
    set(findobj(af, 'type', 'line', 'Tag', 'Median'), 'linewidth',4);
end
ax = gca;
if max(groupx(:)) <= 1
    ax.XTick = unique(groupx(:));
    ax.XTickLabel = strcat(string(round(100*sort(ax.XTick))),'%');%strcat(string([(ax.XTick)*100/unitb]),'%');
    xlabel('Time to storm origin [%]')
    % xlim([0.5,numel(unique(groupx(:)))+0.5]);
else
    % ax.XTickLabel = string(round(sort(unique(groupx(:)))));
    ax.XTick = [0,12,24,36,48,60,72,84,96];
    xlabel('Event duration (h)')
end
ax.YTick = [0,0.2,0.5,0.8,1];
ylabel('val','FontSize',8)
set(gca,'linewidth',1)
grid off
box off
ax = gca;
set(gca,'TickDir','out');
xtickangle(90)

if 0
    colback = copper(20);
    for loc = 1:numel(unique(groupx(:)))
        patch([0.5,1.5,1.5,0.5]+(loc-1),reshape(repmat(ax.YLim,2,1),1,[]),colback(loc,:),'EdgeColor','none');
        hold on;
    end
    alpha(0.1)
end

ylabel(ylabelText);
% box on;

end

function valthis = plotEvalOne(EXPS,evalM,ylab)
% Example:
%     EXPS = {'C1U5-1','C1U10-1','C1P5-Rad1','C1P5-1','C1P5In-Rad1','C1P5In-1'};
%     evalM = 'hit';%'fa','cci','eb'
%     plotEvalOne(EXPS,evalM,ylab)
savePath = ['K:\DATA_FCS\CrossVal\'];
[col, ~, ~, ~, ~] = cptcmap('cw1-002','mapping','scaled','ncol',10);
% col = bone(10);
mark = {'o','s','^','d','p','>','.'};
tag = 1;
h = handle(0);
hName = [];
[HRval,FAval,CSIval,EBval] = deal([]);
for expNo = EXPS
    
    expNo = string(expNo);
    % [version,testConfig,filefolder] = getExpNoInfo(expNo);
    
    for configVec = [0;0;0]%[0 1 0 0; 0 0 1 0; 0 0 0 1]
        
        testConfig.IDW = logical(configVec(1));% false is better (at least for event 90)
        testConfig.uniqueEvent = logical(configVec(2));
        testConfig.includePrevious = logical(configVec(3));
        load([savePath,sprintf('SummarySTATS_%s_IDW%g_uniqueEvent%g_includePrevious%g.mat',...
            expNo,testConfig.IDW,testConfig.uniqueEvent,testConfig.includePrevious)],...
            'STATS','testConfig');
        % STATS(STATS.eventTi>STATS.eventTT*0.8 | STATS.eventTi<STATS.eventTT*0.2,:) = [];
        
        statarray = grpstats(STATS,'eventNo',{'median'});
        eval(['valthis(:,tag) = statarray.median_',evalM,';']);
        % try
        %     eval(['valthis(:,tag) = STATS.',evalM,';']);
        % catch me
        %     eval(['valthis(1:size(STATS.',evalM,',1),tag) = STATS.',evalM,';']);
        %     eval(['valthis(size(STATS.',evalM,',1)+1:end,tag) = NaN;']);
        % end
    end
    tag = tag+1;
end


setFigureProperty('Subplot2');
set(0,'defaultAxesFontSize',10,'defaultAxesFontName','Arial',...
    'defaultAxesTitleFontWeight','Bold',...
    'DefaultLineLineWidth', 2,...
    'defaultTextFontSize',10);
boxplot(valthis(:), reshape(repmat(EXPS(1,:),size(valthis,1),1),1,[]),'colors',col,'Whisker',1000,...
    'symbol','ko');
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(findobj(gcf,'LineStyle','--'),'LineStyle','-')
set(lines, 'linewidth',4);

ylabel('val','FontSize',8)
set(gca,'linewidth',1)
grid off
ax = gca;
ax.XTick = 1:length(EXPS);
ax.XTickLabel = EXPS;
set(gca,'TickDir','out');
xtickangle(90)

if strcmp(evalM,'eb')
    set(ax,'YScale','log');
    ylim([1e-3,1000]);
    ax.YTick = [0.01,0.1,1,10,100,1000];
    ax.YTickLabel = [ax.YTick(1:end-1),Inf];
else
    set(ax,'YScale','linear');
    ylim([-0.1,1.1])
end

% colback = pink(length(EXPS));
% for loc = 1:length(EXPS)
%     patch([0.5,1.5,1.5,0.5]+(loc-1),reshape(repmat(ax.YLim,2,1),1,[]),colback(loc,:),'EdgeColor','none');
%     hold on;
% end

hold on;
plot(ax.XLim,[1,1],'k--','linewidth',1);
alpha(0.1)
ylabel(ylab);
box off

end


function [ETab,TP,FN,FP,TN,evi] = plotEvalMat_Mean(EXPS)
savePath = ['K:\DATA_FCS\CrossVal\'];

hName = [];
% [HRval,FAval,CSIval,EBval] = deal([]);
[TPRval,TNRval,PPRval,NPRval,ACCval] = deal([]);
for expNo = EXPS
    
    expNo = string(expNo);
    % [version,testConfig,filefolder] = getExpNoInfo(expNo);

    for configVec = [0;0;0]%[0 1 0 0; 0 0 1 0; 0 0 0 1]
        
        testConfig.IDW = logical(configVec(1));% false is better (at least for event 90)
        testConfig.uniqueEvent = logical(configVec(2));
        testConfig.includePrevious = logical(configVec(3));
        load([savePath,sprintf('SummarySTATS_%s_IDW%g_uniqueEvent%g_includePrevious%g.mat',...
            expNo,testConfig.IDW,testConfig.uniqueEvent,testConfig.includePrevious)],...
            'STATS','testConfig');
        % STATS(STATS.eventTi>STATS.eventTT*0.8 | STATS.eventTi<STATS.eventTT*0.2,:) = [];
        prc25 = @(x)prctile(x,25);
        prc75 = @(x)prctile(x,75);
        getM = @(x)nanmedian(x);%nanmedian(x);
        
        statarray = grpstats(STATS,'eventNo',{'median','mode','std',@(x)prc25(x),@(x)prc75(x), @(x)std(x)./nanmean(x)});
        ETab = [getM(statarray.Fun4_tpr),getM(statarray.median_tpr),getM(statarray.Fun5_tpr),getM(statarray.Fun6_tpr);...
            getM(statarray.Fun4_tnr),getM(statarray.median_tnr),getM(statarray.Fun5_tnr),getM(statarray.Fun6_tnr);...
            getM(statarray.Fun4_ppr),getM(statarray.median_ppr),getM(statarray.Fun5_ppr),getM(statarray.Fun6_ppr);...
            getM(statarray.Fun4_npr),getM(statarray.median_npr),getM(statarray.Fun5_npr),getM(statarray.Fun6_npr);...
            getM(statarray.Fun4_acc),getM(statarray.median_acc),getM(statarray.Fun5_acc),getM(statarray.Fun6_acc)];
        
        [TP,FN,FP,TN] = deal(NaN(size(STATS,1),1));
        for imn = 1:size(STATS,1)
            [TP(imn),FN(imn),FP(imn),TN(imn)] = SolveIt(STATS.tpr(imn),STATS.tnr(imn),...
                STATS.ppr(imn),STATS.npr(imn),STATS.acc(imn));
        end
        evi = STATS.eventNo;
    end

end

end

function [x1,x2,x3,x4] = SolveIt(y1,y2,y3,y4,y5)
% Y: tpr,tnr,ppr,npr,acc
% X: tp,fn,fp,tn

A = [1-y1,-y1,0,0;
    0,0,-y2,1-y2;
    1-y3,0,-y3,0;
    0,-y4,0,1-y4;
    1,0,0,1];
b = [0,0,0,0,140*y5]';
X = A\b;
x1 = X(1);x2 = X(2);x3=X(3);x4=X(4);
end

function ETab = plotEvalMat_50Prc(EXPS)
savePath = ['K:\DATA_FCS\CrossVal\'];
col = bone(6);
mark = {'o','s','^','d','p'};
tag = 1;
h = handle(0);
hName = [];
% [HRval,FAval,CSIval,EBval] = deal([]);
[TPRval,TNRval,PPRval,NPRval,ACCval] = deal([]);
for expNo = EXPS
    
    expNo = string(expNo);
    % [version,testConfig,filefolder] = getExpNoInfo(expNo);
    
    col = copper(6);%
    %     (tag==1) * copper(6) + ...
    %         (tag==2) * gray(6)+...
    %         (tag==3) * bone(6)+...
    %         (tag==4) * pink(6)+...
    %         (tag==5) * hsv(6)+...
    %         (tag==6) * hot(6);
    xloc = (tag==1) * ((1:5)-0.2) +...
        (tag==2) * (1:5) +...
        (tag==3) * ((1:5)+0.2) +...
        (tag==4) * ((1:5)+0.4) +...
        (tag==5) * ((1:5)+0.6) +...
        (tag==6) * ((1:5)+0.8);
    
    for configVec = [0;0;0]%[0 1 0 0; 0 0 1 0; 0 0 0 1]
        
        testConfig.IDW = logical(configVec(1));% false is better (at least for event 90)
        testConfig.uniqueEvent = logical(configVec(2));
        testConfig.includePrevious = logical(configVec(3));
        load([savePath,sprintf('SummarySTATS_%s_IDW%g_uniqueEvent%g_includePrevious%g.mat',...
            expNo,testConfig.IDW,testConfig.uniqueEvent,testConfig.includePrevious)],...
            'STATS','testConfig');
        % STATS(STATS.eventTi>STATS.eventTT*0.8 | STATS.eventTi<STATS.eventTT*0.2,:) = [];
        prc25 = @(x)prctile(x,25);
        prc75 = @(x)prctile(x,75);
        getM = @(x)nanmedian(x);
        
        statarray = grpstats(STATS,'eventNo',{'median','mode','std',@(x)prc25(x),@(x)prc75(x), @(x)std(x)./nanmean(x)});
        % {'mean','median','mode','std','min','max','skewness','kurtosis'}
        
        % ETab = [getM(statarray.Fun4_hit),getM(statarray.median_hit),getM(statarray.Fun5_hit),getM(statarray.Fun6_hit);...
        %     getM(statarray.Fun4_fa),getM(statarray.median_fa),getM(statarray.Fun5_fa),getM(statarray.Fun6_fa);...
        %     getM(statarray.Fun4_cci),getM(statarray.median_cci),getM(statarray.Fun5_cci),getM(statarray.Fun6_cci);...
        %     getM(statarray.Fun4_eb),getM(statarray.median_eb),getM(statarray.Fun5_eb),getM(statarray.Fun6_eb)];
        ETab = [getM(statarray.Fun4_tpr),getM(statarray.median_tpr),getM(statarray.Fun5_tpr),getM(statarray.Fun6_tpr);...
            getM(statarray.Fun4_tnr),getM(statarray.median_tnr),getM(statarray.Fun5_tnr),getM(statarray.Fun6_tnr);...
            getM(statarray.Fun4_ppr),getM(statarray.median_ppr),getM(statarray.Fun5_ppr),getM(statarray.Fun6_ppr);...
            getM(statarray.Fun4_npr),getM(statarray.median_npr),getM(statarray.Fun5_npr),getM(statarray.Fun6_npr);...
            getM(statarray.Fun4_acc),getM(statarray.median_acc),getM(statarray.Fun5_acc),getM(statarray.Fun6_acc)];
        
        % HRval(tag,:) = ETab(1,:);
        % FAval(tag,:) = ETab(2,:);
        % CSIval(tag,:) = ETab(3,:);
        % EBval(tag,:) = ETab(4,:);
        TPR(tag,:) = ETab(1,:);
        TNR(tag,:) = ETab(2,:);
        PPR(tag,:) = ETab(3,:);
        NPR(tag,:) = ETab(4,:);
        ACC(tag,:) = ETab(5,:);
        
        setFigureProperty();
        % ETab(ETab>=5) = 5;
        h1 = plot(xloc,ETab(:,2),mark{tag},'color',col(tag,:),...
            'markerfacecolor',[0.5 0.5 0.5],'markersize',isequal(configVec,[0;0;0])*12+~isequal(configVec,[0;0;0])*2);
        
        hold on
    end
    h(tag) = h1(1);
    hName{tag} = expNo;
    tag = tag+1;
end
%%
set(gca,'linewidth',1)
grid off

ax = gca;xlim([0.5 5.5]);ylim([0.1,1.1])
ax.XTick = 1:5;

ax.XTickLabel = {'TPR','TNR','PPR','NPR','ACC'};
ax.YTick = 0:0.2:1;
ax.YTickLabel = 0:0.2:1;
% ax.XTickLabel = {'HR','FA','CSI','EB'};
% ax.YTick = [0.2:0.2:0.8,1:0.5:4];
% ax.YTickLabel = [0.2:0.2:0.8,1:0.5:4];

patch([0.5,1.5,1.5,0.5],reshape(repmat(ax.YLim,2,1),1,[]),'r','EdgeColor','none');
patch([0.5,1.5,1.5,0.5]+1,reshape(repmat(ax.YLim,2,1),1,[]),'y','EdgeColor','none');
patch([0.5,1.5,1.5,0.5]+2,reshape(repmat(ax.YLim,2,1),1,[]),'b','EdgeColor','none');
patch([0.5,1.5,1.5,0.5]+3,reshape(repmat(ax.YLim,2,1),1,[]),'g','EdgeColor','none');
patch([0.5,1.5,1.5,0.5]+4,reshape(repmat(ax.YLim,2,1),1,[]),'m','EdgeColor','none');

plot(ax.XLim,[1,1],'k--');
alpha(0.1)
legend(h,hName,'location','southwest');

set(ax,'YScale','linear');
ylabel('50 Prctile')

end

function formatHist2(X,group,ylabelText,edges)

EDGES = [edges,Inf];%[0,0.5,1,2,Inf];%[0:0.1:2,Inf];
eventNo = 1:numel(unique(group));
GMAP = NaN(numel(EDGES)-1,length(eventNo));

for evi = eventNo
    thisX = X(group == evi);
    [H,~] =  histcounts(thisX,EDGES);
    GMAP(:,evi) = H./sum(H);
end
GMAP(GMAP==0)=NaN;
X = 1:length(EDGES(1:end-1));
Xtrue = EDGES(1:end-1);
Y = eventNo;

Y0 = Y;% 1:length(Y);
pcolor(Y0,Xtrue,GMAP);shading flat
alpha(0.7)
cptcmap('diff_4_bias', 'mapping', 'scaled','ncol',9);
set(gca,'clim',[0,1]);
box on

set(gca,'XTick',Y0)
set(gca,'XTickLabel',Y);
ylabel(ylabelText)
xlabel('eventNo')
colorbar
ax = gca;
ax.XTickLabelRotation = 90;
set(gca,'XTick',[])

end

function saveSummaryStats(expNo,savePath)
for configVec = [0 1 0 0;
                0 0 1 0;
                0 0 0 1]%[0;0;0]
    expNo = string(expNo);
    
    [version,testConfig,filefolder] = getExpNoInfo(expNo);
    testConfig.IDW = logical(configVec(1));% false is better (at least for event 90)
    testConfig.uniqueEvent = logical(configVec(2));
    testConfig.includePrevious = logical(configVec(3));
    
    tstart = tic;
    tend = [];
    for eventNo4val = 1:157
        
        if eventNo4val == 1
            
            load([filefolder,'STATS_',sprintf('%s_ev%03d.mat',version,eventNo4val)],...
                'STATS','version','FOREOUTPUT'); %#ok<NASGU>
            
            [STATS] = evaluateFloodPrediction_FCSVersion(FOREOUTPUT.FLO_test, FOREOUTPUT.FLO_pred,...
                FOREOUTPUT.E, FOREOUTPUT.N, FOREOUTPUT.analogDistL1, testConfig);
            % [STATS] = evaluateFloodPrediction(FOREOUTPUT.FLO_test, FOREOUTPUT.FLO_pred,...
            %     FOREOUTPUT.E, FOREOUTPUT.N, FOREOUTPUT.analogDistL1, testConfig);

            STATS = rmfield( STATS, 'eventNo_pred');
            STATS = structfun(@(x)reshape(x,[],1), STATS, 'UniformOutput', false);
            STATS = struct2table(STATS);
            
        else
            
            A = load([filefolder,'STATS_',sprintf('%s_ev%03d.mat',version,eventNo4val)],...
                'STATS','version','FOREOUTPUT');
            [A.STATS] = evaluateFloodPrediction_FCSVersion(A.FOREOUTPUT.FLO_test, A.FOREOUTPUT.FLO_pred,...
                A.FOREOUTPUT.E, A.FOREOUTPUT.N, A.FOREOUTPUT.analogDistL1, testConfig);
            % [A.STATS] = evaluateFloodPrediction(A.FOREOUTPUT.FLO_test, A.FOREOUTPUT.FLO_pred,...
            %     A.FOREOUTPUT.E, A.FOREOUTPUT.N, A.FOREOUTPUT.analogDistL1, testConfig);
            A.STATS = rmfield( A.STATS, 'eventNo_pred');
            A.STATS = structfun(@(x)reshape(x,[],1), A.STATS, 'UniformOutput', false);
            A.STATS = struct2table(A.STATS);
            
            STATS = [STATS;A.STATS];
        end
        % plot. can be deleted. %
        tend = [tend,toc(tstart)];
        semilogy(tend,'r-');set(gca,'linewidth',1)
        title(['ExpNo.',expNo]);xlabel('EventNo');ylabel('time');
        grid minor
        drawnow
        % plot. can be deleted. %
    end
    
    clear eventNo4val A
    
    [startTime,endTime,eventNos] = getHistEventTime(testConfig.rainfallSource);
    eventTime = datetime(datevec(startTime));
    
% load([savePath,filesep,sprintf('SummarySTATS_%s_IDW%g_uniqueEvent%g_includePrevious%g.mat',...
%         expNo,testConfig.IDW,testConfig.uniqueEvent,testConfig.includePrevious)],...
%         'STATS','testConfig');
%     STATS.Properties.VariableNames{7}='ppr';
    save([savePath,filesep,sprintf('SummarySTATS_%s_IDW%g_uniqueEvent%g_includePrevious%g.mat',...
        expNo,testConfig.IDW,testConfig.uniqueEvent,testConfig.includePrevious)],...
        'STATS','testConfig');
end
end
