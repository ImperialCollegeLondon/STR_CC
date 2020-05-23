% This script is to plot IDAF based on annual maxima ananlysis
%
% Data: 1980-2000 UKCP CPM rainfall or 2060-2080 UKCP CPM rainfall
%
%
%
% ref: ...2002
% based on C. De Michele 2011 func.6 in Journal of Hydrology
%
% @ Yuting CHEN
% yuting.chen17@imperial.ac.uk
%


clear;clc
close all


futPeriod = '2060-2080';
pasPeriod = '1980-2000';
dataInfo = getDataInfo('SCO',pasPeriod);
rt = exp(linspace(log(1.2),log(40),80));
area = 2.^(unique(round(5*log2(1:8000)))/5);
figForm = struct; % format of output figure.
[figForm.rt,figForm.area] = meshgrid(rt, area);
itag = 1; Zs = [];
REGIONS = {'SCO','WAL','EUK'};

ha = tight_subplot(3,1,[.0 .0],[.15 .05],[.30 .05]);

for regionName = REGIONS
    
    regionName = regionName{1};
    for dataInfo = [getDataInfo(regionName,pasPeriod),getDataInfo(regionName,futPeriod)]
        Imap = [];
        try
        for ensNo = dataInfo.ensNo(1:end) % dataInfo.ensNo(:)'
            ensNo = ensNo{1};
            load(sprintf('%s%sIDAFs_%s_ensNo%s_%s.mat',dataInfo.savePath,filesep,...
                dataInfo.region.Name,ensNo,dataInfo.Years),'T','TIMES');
            pl = 0;% not plot
            parmhat_gev = grpstats(T.am,T.area,{@(x)fitAM(x,pl)});
            
            AREA = grpstats(T.area,T.area,{'mean'});
            RT = 1.1:0.1:40;
            I0 = getI(RT,parmhat_gev);
            [rt0,area0] = meshgrid(RT,AREA);
            Imap = cat(3,Imap,I0);

        end
        catch me
            disp(ensNo);
        end
        
        if strcmp(dataInfo.Years,pasPeriod)
            pasImap = Imap;
        elseif strcmp(dataInfo.Years,futPeriod)
            futImap = Imap;
        end
        
%{
        I = nanmedian(Imap,3);
        I = griddata(rt0,area0,I,figForm.rt,figForm.area);
        figure;setFigureProperty('Paper')
        hold on;
        pcolor(rt,log2(area),I);shading flat;
        [h,c] = contour(rt,log2(area),I,[20:20:100],'w-','ShowText','on');
        clabel(h,c,'Color',[0.9 0.9 0.9])
        c.LineWidth = 1;
        cptcmap('cw1-002','mapping','scaled','ncol',100);
        set(gca,'YScale','linear')
        ax = gca;
        ax.YLim = log2([5,8000]);
        ax.YTick = log2([5,10,20,50,100,200,500,1000,2000,4000,8000]);
        ax.YTickLabels = 2.^ax.YTick;
        ax.TickDir = 'out';
        caxis([0,150]);
        colorbar
        title(sprintf('%s-%s',regionName,dataInfo.Years))
        xlabel('Return Period (yr)')
        ylabel('Catchment Area (km^2)')
        savePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\IDAF\';
        saveName = sprintf('%s-IDAF-%s-%s-Dur%02dh',regionName,dataInfo.Years,'allEns',1);
        savePlot([savePath,filesep,saveName],'units','points','XYWH',[150,0,300,230],'needreply','N','onlyPng',true);
%}        
    end
    
    %%
    
    prcval = 50;
    Z = (futImap-pasImap)./pasImap;% Z = nanmedian(futImap-pasImap,3);
    Zs{itag} = Z;
    Z = 100*prctile(Z,prcval,3);%nanmedian(Z,3);
    
    sigZ = testSig(permute(futImap,[3,1,2]),permute(pasImap,[3,1,2]));
    Z = griddata(rt0,area0,Z,figForm.rt,figForm.area);
    sigZ = griddata(rt0,area0,sigZ,figForm.rt,figForm.area,'nearest');
    
    
    axes(ha(itag))
    setFigureProperty('Paper')
    set(0,'defaultAxesFontSize',10);
    % hold on;
    pcolor(rt,log2(area),Z);shading flat;hold on;
    % alpha(0.95)%'diff_darkBlue_darkRed'%BlueDarkRed18
    cptcmap('RdYlBu_11','mapping','scaled','ncol',11,'flip',false);% ncol 20
    % [~,c] = contour(rt,log2(area),Z,[-30,-20,-10,10,15,20,25,30],'k-','ShowText','on');
    [h,c] = contour(rt,log2(area),Z,[-40,-20,-10,10,20,30],'w-','ShowText','on');
    clabel(h,c,'Color',[0.95 0.95 0.95])
    plot(figForm.rt(mod(figForm.rt,2)>=0 & sigZ==0),...
        log2(figForm.area(mod(figForm.rt,2)>=0 & sigZ==0)),'.','color',[0.2 0.2 0.2 0.5],...
        'markersize',2);
    c.LineWidth = 1;
    ax = gca;
    ax.YScale = 'linear';
    ax.XScale = 'log';
    ax.YLim = log2([5,8000]);
    ax.XLim = [1.2,39.9];
    % ax.YTick = log2([5,10,20,50,100,200,500,1000,2000,4000,8000]);
    ax.YTick = log2([20,50,100,200,500,1000,3000,8000]);
    ax.YTickLabels = 2.^ax.YTick;
    ax.XTick = [2,5,10,20,40];
    ax.TickDir = 'out';
    caxis([-35,35]);
    box off
    % c = colorbar;
    % c.Ticks = [-50,-30,-10,0,10,30,50];
    % c.TickLabels = strcat(c.TickLabels,'%');
    % title(sprintf('(future-past)/past'))
    xlabel('RT (yr)')
    if itag == 2
    ylabel('Area (km^2)')
    end
    savePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\IDAF\';
    
    itag = itag+1;
end

c = colorbar('SouthOutside');
% c.Ticks = [-50,-30,-15,0,15,30,50];
c.Ticks = [-30,-15,0,15,30];
c.TickLabels = strcat(c.TickLabels,'%');

% saveName = sprintf('%s-diff_IDAF-%s-Dur%02dh_percentile%02d',regionName,'allEns',1,prcval);
saveName = sprintf('diff_IDAF-%s-Dur%02dh_percentile%02d','allEns',1,prcval);
savePlot([savePath,filesep,saveName],'units','centimeters',...
    'XYWH',[5,0,5,12],'needreply','Y','onlyPng',false);
close all
%%
% easier-understanding plot
figure;
setFigureProperty('Subplot2');
set(0,'defaultAxesFontSize',10);
hold on;
h = handle(0);
ax = gca;% xlim([2,(8000)])
ax.XLim = log2([5.3,8000]);
ax.YLim = [-6,41];
ax.XTick = log2([5,10,20,50,100,200,500,1000,2000,4000,8000]);
ax.XTickLabels = 2.^ax.XTick;
% patch(log2([20,500,500,20,20]),[ax.YLim([1,1,2,2,1])],[1 0 0],...
%     'edgecolor','none');

col = {[66 146 199]/255,[240  60  43]/255,[0.3 0.3 0.3]};
for i = [1,3,2] %regionName = {'SCO','WAL','EUK'}
    % subplot(1,3,i);hold on;
    Zp50 = 100*prctile(Zs{i},50,3);%
    Zmean = 100*nanmean(Zs{i},3);
    Z0 = Zp50.*(repmat(200./RT,size(Zp50,1),1))/nansum(200./RT);
    h(i) = plot(log2(AREA),nansum(Z0,2),'-','Color',col{i},'linewidth',1);
    Z0_mean = Zmean.*(repmat(200./RT,size(Zmean,1),1))/nansum(200./RT);
    % plot(log2(AREA),nansum(Z0_mean,2),'o','Color',col{i},'linewidth',1);
    % close all

    Zp20 = 100*prctile(Zs{i},20,3);%
    Z10 = Zp20.*(repmat(200./RT,size(Zp20,1),1))/nansum(200./RT);
    Zp80 = 100*prctile(Zs{i},80,3);%
    Z90 = Zp80.*(repmat(200./RT,size(Zp80,1),1))/nansum(200./RT);
    fill([log2(AREA);flip(log2(AREA))],[nansum(Z90,2);flip(nansum(Z10,2))],h(i).Color,'LineStyle','none');
    alpha(0.12/i)
%     for perval = linspace(1,100,12)
%     Zp_1 = 100*prctile(Zs{i},perval,3);%
%     Z_1 = Zp_1.*(repmat(200./RT,size(Zp_1,1),1))/nansum(200./RT);
%     plot(log2(AREA),nansum(Z_1,2),'Color',h(i).Color,'linewidth',0.5)
%     end
    % close all
end
% title('Change of I(A,1h,RT)')

xlabel('Catchment Area /Km^2');
ylabel('Relative Change(%)');
set(gca,'TickDir','Out')
xtickangle(45)
% box on;

legend(h,getPrintName(REGIONS))

savePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\IDAF\';
saveName = sprintf('ExtremeChange(CatchmentArea)');
savePlot([savePath,filesep,saveName],'units','centimeters',...
    'XYWH',[5,0,12,12.3],'needreply','Y','onlyPng',false);

%% Auxillary function
function [isSigDif] = testSig(Z1,Z2)
% Z: [loc1,loc2,sampleSize]

isSigDif = squeeze(ttest(Z1,Z2,'alpha',0.1));

end
function [pname] = getPrintName(REGIONS)
if iscell(REGIONS)
    pname = [];
    fprintf('Make sure <REGIONS> is in this order: SCO,WAL,EUK!\n');
    pname = {'NUK','SWUK','SEUK'};
end
end
function [I] = getI(RT,parmhat)
I = [];
P = 1-(1./RT);
for leni = 1:size(parmhat,1)
    I(leni,:) = gevinv(P,parmhat(leni,1),parmhat(leni,2),parmhat(leni,3));
end
end
function [parmhat] = fitAM(a,pl)
% get stats for a series of a;

[parmhat,parmci] = gevfit(a);

x = (.5:1:(numel(a)-.5))'./ numel(a);
r = 1./(1-x);
T = 1.05:0.05:40;
P = 1-(1./T);

X = gevinv(P,parmhat(1),parmhat(2),parmhat(3));

if pl
    colorI = 'r';
    hobs = plot(r,sort(a),'ko','markerfacecolor',colorI,'markersize',5); % cdf of y
    hold on;
    hfit = plot(T,X,'color',colorI);
    xlim([0,40])
    pause(0.02);
    hold off;
end
parmhat = parmhat';
end


function dataInfo = getDataInfo(regionName,regionYear)

dataInfo = struct;
dataInfo.Years = regionYear;
dataInfo.region = getfield(REGIONS_info(),regionName);
dataInfo.ensNo = getEnsNos();
if strcmp(regionYear,'2060-2080')
    dataInfo.fileGetPath = ['D:/UKCP18_Future_2060_2080/',dataInfo.region.Name];
    dataInfo.savePath = 'D:/UKCP18/IDAF';
elseif strcmp(regionYear,'1980-2000')
    dataInfo.fileGetPath = ['D:/UKCP18/',dataInfo.region.Name];
    dataInfo.savePath = 'D:/UKCP18/IDAF';
end
dataInfo.figPath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP';
dataInfo.season = 'JJA';
dataInfo.scaleF = 32;
dataInfo.gridReso = 2.2;

idaf = struct;
idaf.duration =  [1];%[ 1, 3, 6,12,24];
idaf.method = 'annualMaximum'; % reference: C. De Michele 2011 in JH
idaf.area = [1:100].^2*(dataInfo.gridReso).^2;
dataInfo.idaf = idaf;

end

