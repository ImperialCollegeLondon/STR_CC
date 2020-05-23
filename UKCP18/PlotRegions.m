% This file stores all scripts for plotting the regional patterns.
% PlotRegions
%
% @ Yuting Chen
% Imperial College London
% Update: 2019.12.06
% 
% Notice here the original matrix order of CPMs is [E,N,T]
%                                          radar is [N,E,T];
% So conversion (permute) was undertaken in func loadRainNimrod() to 
%    change the output order of CPMs into [E,N,T];
                                         
                                         
colm = pink(8);
colm(1:2,:) = [];
colOBS = copper(8);
colOBS(1:2,:) = [];
lineType = 'x';%'-';
lineTypeObs = '^';% 'o-';
linewidth = 1;
linewidthObs = 2;
Name = {'London','SWestuk','Westuk','Scotland'};%{'London','SWestuk','Westuk','Scotland'};
dataSP = 'H:\DATA_CLIMATE\UKCP18\';

REGIONS = REGIONS_info();

%% Figure .: Hourly rainfall for given percentile threshold
% This figure follows the processes in Kendon 2012 Fig.3
% in order to compare the difference of cdf tail between radar and 2.2km
% cpm

figure;
setFigureProperty('Paper');
XYWH = [50,-50,700,500];
set(gcf,'units','points','position',XYWH);
ha = tight_subplot(2, 2,[.1 .1],[.1 .1],[.1 .1]);

region = REGIONS.Westuk; %London;%SWestuk; %Scotland;

% for cpm
Period = '1980-2000';
[dataplot] = cell(1,4);
colm = pink(20);
for thissp = 2%1:4
    
    season = thissp;
    
    axes(ha(thissp));
    
    [E,N,RainEnsembles,scaleF,region] = loadRainEnsembles(region,getMons(season),Period);
    [hsim] = plotPrcRain(E,N,RainEnsembles,scaleF,region,colm);
    
    fprintf('------thissp_%1d------\n',thissp)
    
    % xlim([0.85,0.9999]);
    ylim([0.1,20]);
    set(gca,'YScale','log')
    
end

% for radar
[dataplot] = cell(1,4);

for thissp = 2%1:4
    
    axes(ha(thissp));
    
    season = thissp;

    [E,N,RainNimrod,scaleF,region] = loadRainNimrod(region,getMons(season));
    [hobs] = plotPrcRain(E,N,{RainNimrod},scaleF,region,[1 0 0]);
    
    fprintf('------thissp_%1d------\n',thissp)
end

for thissp = 2%1:4
    season = thissp;
    axes(ha(thissp));
    % xlim([0.85,0.9999]);
    ylim([0.02,20]);
    grid minor;
    set(gca,'YTick', yticks)
    set(gca,'TickLength',[0.03, 0.018])
    legend(gca,[hsim,hobs],{'2.2km','radar'},'Location','NorthWest')
    title(getSeasonName(season))
    % text(99.99,20,getSeasonName(season),'verticalalignment','top',...
    %     'horizontalalignment','right');
    xlabel('percentile');
    ylabel('mm/h')
    xtickangle(90);
    ax = gca;
    box off
    ax.YTick = [0.1,0.2,0.5,1,5,10,20,40];
    ax.YTickLabel = [0.1,0.2,0.5,1,5,10,20,40];
    set(gca,'TickDir','out');
end

format_xylabel(ha,2,2)

savePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP';
saveName = sprintf('%s_Kendon2012RainPrctile',region.Name);
savePlot([savePath,filesep,saveName],'XYWH',[150,0,700,500],'needreply','N','onlyPng',true);
%% Figure 0: Spatial Spectrum Analysis
region = REGIONS.London;%SWestuk; %Westuk; %Scotland;
pltag = 1;
reg = {REGIONS.London,REGIONS.SWestuk,REGIONS.Westuk,REGIONS.Scotland};
% {REGIONS.London,REGIONS.Westuk,REGIONS.Scotland};
% 
colm = pink(20);
figure;
setFigureProperty('Paper');
set(0,'defaultAxesFontSize',10,'defaultTextFontSize',10);
set(gcf,'units','centimeters','position',[5,0,18,10]);
ha = tight_subplot(2,length(reg),[.01 .01],[.1 .05],[.1 .01]);

for seasonInd = [2,4]% 1:4
    for regind = 1:length(reg)
        
        region = reg{regind};
        
        axes(ha(pltag));
        mons = getMons(seasonInd);
        PeriodName = getSeasonName(seasonInd);
        
        ggg(region,mons,PeriodName,colm,colOBS,Period);
        
        pltag = pltag+1;

        drawnow
        
    end
end
format_xylabel(ha,2,length(reg))
filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\';
filename = [filePath,filesep,'SpatialSpectrumAnalysis'];
savePlot(filename,'units','centimeters','XYWH',[5,0,18,10],'needreply','Y');

%% FIGURE; SCATTER plot

%% Figure 1: Plot the saved dataplot

[h,hobs] = deal([]);

figure;
XYWH = [50,-50,500,500];
set(gcf,'units','points','position',XYWH);
ha = tight_subplot(4,4,[.00 .00],[.1 .01],[.1 .1]);

    
logsf = 8;

for thisstation = 1:length(Name)
    CPM = load(sprintf('%sWAR_Area_%s_4Season.mat',dataSP,Name{thisstation}),'dataplot');  
    RAD = load(sprintf('%sObs_WAR_Area_%s_4Season.mat',dataSP,Name{thisstation}),'dataplot');
    
    tsp=1;
    for thisseason = [2 3 4 1]    
        thisplot = (thisstation-1)*4+tsp;%thisseason;
        tsp=tsp+1;
        axes(ha(thisplot));
        
        
        indThr = 1:24;
        y = CPM.dataplot{thisseason}.y1(:,indThr);
        x = repmat(RAD.dataplot{thisseason}.y1(indThr),[12,1]);
        plot(x'*100,y'*100,'k:','color',0.5*ones(1,3),'linewidth',1);
        hold on;
        
        indThr = [1,2,6,8,11,16];%1:24;
        y = CPM.dataplot{thisseason}.y1(:,indThr);
        x = repmat(RAD.dataplot{thisseason}.y1(indThr),[12,1]);
        c = repmat(flip(hsv(numel(indThr))),[1,1,12]);
        c = reshape(permute(c,[3,1,2]),[],3);
        siz = repmat((indThr)*3+20,[12,1]);
        
        scatter(logsf*reshape(x*100,[],1),logsf*reshape(y*100,[],1),siz(:),c,'^','filled','linewidth',1,...
            'markeredgecolor','k')
        
        hold on;
        plot([0,0.000002,15]*logsf,[0,0.000002,15]*logsf,'k-','linewidth',0.5)
        xlabel('RAD');
        ylabel('CPM2.2')
        set(gca,'linewidth',1)
        box on
        text(0.0002*logsf,15*logsf,sprintf('%s-%s',Name{thisstation},getSeasonName(thisseason)),'fontsize',12,...
            'horizontalalignment','left',...
            'verticalalignment','top');
        %{
        [h(thisseason),~] = plot_WetArea(CPM.dataplot{thisseason}.Thres,...
            CPM.dataplot{thisseason}.y1,CPM.dataplot{thisseason}.y2,...
            colm(thisseason,:),lineType,linewidth);
        [hobs(thisseason),~] = plot_WetArea(RAD.dataplot{thisseason}.Thres,...
            RAD.dataplot{thisseason}.y1,RAD.dataplot{thisseason}.y2,...
            colOBS(thisseason,:),lineTypeObs,linewidthObs);
        %}
 
       
        if 1%thisstation == length(Name)
            
            grid minor
            set(gca,'XScale','log')
            set(gca,'YScale','log')
            ylim([0.0002*logsf,15*logsf]);
            xlim([0.0002*logsf,15*logsf])
            ytickformat(gca, '%.2f%%');
            xtickformat(gca, '%.2f%%');
            ytickformat('percentage')
            yticks([1e-5 1e-4 1e-3 1e-2 1e-1 1 10]*logsf)
            yticklabels({'1e-5%','1e-4%','1e-3%','0.01%','1%','10%'})
            xtickformat('percentage')
            xticks([1e-5 1e-4 1e-3 1e-2 1e-1 1 10]*logsf)
            xticklabels({'1e-5%','1e-4%','1e-3%','0.01%','1%','10%'})
            
            % text(0,0.999,getSeasonName(thisseason),'fontsize',12);
            ax = gca;
            ax.FontSize = 12;
            ax.XTickLabelRotation = 90;
            % legend boxoff
            
            % cptcmap('diff_darkBlue_darkRed', 'mapping','scaled');
            
        end
    end
end
XYWH = [50,-50,500,500];
set(gcf,'units','points','position',XYWH);

% c = colorbar('location','Manual', 'position', [0.85 0.1 0.02 0.30]);
% c.Ruler.TickLabelFormat='%g%%';

format_xylabel(ha,4,4)

% plot_WetArea(dataplot.Thres,dataplot.y1,dataplot.y2,colm(season,:));

filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\';
filename = [filePath,filesep,'WAR_rad_cpm_4regions_normalAxis'];



%% Figure 1: Plot the saved dataplot

[h,hobs] = deal([]);

figure;
XYWH = [50,-50,500,500];
set(gcf,'units','points','position',XYWH);
ha = tight_subplot(2,2,[.05 .05],[.1 .01],[.1 .1]);
cmap = pink(length(Name)*2);

logsf = 1;%5;
for thisseason = [2 3 4 1]
    tsp=1;
    thisplot = thisseason;% different location in same figure;%(thisstation-1)*4+tsp;%thisseason;
    tsp=tsp+1;
    axes(ha(thisplot));
    for thisstation = 1:length(Name)
        CPM = load(sprintf('%sWAR_Area_%s_4Season.mat',dataSP,Name{thisstation}),'dataplot');
        RAD = load(sprintf('%sObs_WAR_Area_%s_4Season.mat',dataSP,Name{thisstation}),'dataplot');

        indThr = 1:24;
        y = CPM.dataplot{thisseason}.y1(:,indThr);
        x = repmat(RAD.dataplot{thisseason}.y1(indThr),[12,1]);
        % plot(logsf*x'*100,logsf*y'*100,'k:','color',cmap(thisstation,:));%,0.5*ones(1,3),'linewidth',1);
        hold on;
        
        indThr = 1:11;%[1,2,6,8,11];%1:24;
        y = CPM.dataplot{thisseason}.y1(:,indThr);
        x = repmat(RAD.dataplot{thisseason}.y1(indThr),[12,1]);
        c = flip(hsv(numel(indThr)));
        c = reshape(permute(c,[3,1,2]),[],3);
        siz = (indThr)*3+20; % repmat((indThr)*3+20,[12,1]);
        
        scatter(logsf*reshape(nanmean(x,1)*100,[],1),logsf*reshape(nanmean(y,1)*100,[],1),siz(:),c,'o','filled','linewidth',0.5,...
            'markeredgecolor','k')
        
        hold on;
        % plot([0,0.0004,15]*logsf,[0,0.0004,15]*logsf,'k-')
        xlabel('RAD');
        ylabel('CPM2.2')
        set(gca,'linewidth',1)
        box on
       
       
        if thisstation == length(Name)
            XLIM = xlim;
            YLIM = ylim;
            text(0.0004*logsf,15*logsf,sprintf('%s',...%Name{thisstation},...
            getSeasonName(thisseason)),'fontsize',12,...
            'horizontalalignment','left',...
            'verticalalignment','top');
            grid minor
            % set(gca,'XScale','log')
            % set(gca,'YScale','log')
            ylim([min([XLIM,YLIM])*logsf,max([XLIM,YLIM])*logsf]);
            xlim([min([XLIM,YLIM])*logsf,max([XLIM,YLIM])*logsf])
            % ytickformat(gca, '%.2f%%');
            % xtickformat(gca, '%.2f%%');
            % ytickformat('percentage')
            % yticks([1e-2 1e-1 1 10]*logsf)
            % yticklabels({'0.01%','1%','10%'})
            % xtickformat('percentage')
            % xticks([1e-2 1e-1 1 10]*logsf)
            % xticklabels({'0.01%','1%','10%'})
            
            % text(0,0.999,getSeasonName(thisseason),'fontsize',12);
            ax = gca;
            ax.FontSize = 12;
            ax.XTickLabelRotation = 90;
            % legend boxoff
            
            % cptcmap('diff_darkBlue_darkRed', 'mapping','scaled');
            
        end
        
    end
    
end
XYWH = [50,-50,500,500];
set(gcf,'units','points','position',XYWH);

% c = colorbar('location','Manual', 'position', [0.85 0.1 0.02 0.30]);
% c.Ruler.TickLabelFormat='%g%%';

format_xylabel(ha,2,2)

% plot_WetArea(dataplot.Thres,dataplot.y1,dataplot.y2,colm(season,:));

filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\';
filename = [filePath,filesep,'WAR_rad_cpm_4regions_normalAxis'];
savePlot(filename,'XYWH',[150,0,600,600],'needreply','Y');

%% Figure 1: Plot the saved dataplot

[h,hobs] = deal([]);

figure;
XYWH = [50,-50,500,500];
set(gcf,'units','points','position',XYWH);
ha = tight_subplot(4,4,[.00 .00],[.1 .01],[.1 .01]);

for thisstation = 1:length(Name)
    CPM = load(sprintf('%sWAR_Area_%s_4Season.mat',dataSP,Name{thisstation}),'dataplot');  
    for thisseason = 1:4
        thisplot = (thisstation-1)*4+thisseason;%thisseason;
        axes(ha(thisplot));
        [h(thisseason),~] = plot_WetArea(CPM.dataplot{thisseason}.Thres,...
            CPM.dataplot{thisseason}.y1,CPM.dataplot{thisseason}.y2,...
            colm(thisseason,:),lineType,linewidth);
    end
end      
for thisstation = 1:length(Name)
    RAD = load(sprintf('%sObs_WAR_Area_%s_4Season.mat',dataSP,Name{thisstation}),'dataplot');
    for thisseason = 1:4
        thisplot = (thisstation-1)*4+thisseason;%thisseason;
        axes(ha(thisplot));
        [hobs(thisseason),~] = plot_WetArea(RAD.dataplot{thisseason}.Thres,...
            RAD.dataplot{thisseason}.y1,RAD.dataplot{thisseason}.y2,...
            colOBS(thisseason,:),lineTypeObs,linewidthObs);
        if 1%thisstation == length(Name)
            % ylim([0.995,1]);
            xlim([5,35])
            grid minor
            legend(gca,[h(thisseason),hobs(thisseason)],...
                [strcat(Name(thisstation),'-UKCP'),...
                strcat(Name(thisstation),'-Radar')],'Location','SouthEast');
            %legend(gca,[h,hobs],[strcat(Name,'-UKCP'),strcat(Name,'-Radar')],'Location','SouthEast');
            text(15,0.999,getSeasonName(thisseason),'fontsize',12);
            set(gca,'XScale','log')
            legend boxoff
        end
    end
end
XYWH = [50,-50,500,500];
set(gcf,'units','points','position',XYWH);

format_xylabel(ha,4,4)

% plot_WetArea(dataplot.Thres,dataplot.y1,dataplot.y2,colm(season,:));



%% Figure 2: Plot the saved dataplot
figure;
ha = tight_subplot(4,4,[.01 .01],[.1 .01],[.1 .01]);
[h,hobs] = deal([]);
for thisstation = 1:length(Name)
    CPM = load(sprintf('%sIMF_Area_%s_4Season.mat',dataSP,Name{thisstation}),'dataplot');
    RAD = load(sprintf('%sObs_IMF_Area_%s_4Season.mat',dataSP,Name{thisstation}),'dataplot');
    for thisseason = 1:4
        thisplot = (thisstation-1)*4+thisseason;%thisseason;
        axes(ha(thisplot));
        [~,h(thisseason)] = plot_IMF(CPM.dataplot{thisseason}.areaAgg,...
            CPM.dataplot{thisseason}.y1,CPM.dataplot{thisseason}.y2,...
            colm(thisseason,:),linewidth);
        
        [~,hobs(thisseason)] = plot_IMF(RAD.dataplot{thisseason}.areaAgg,...
            RAD.dataplot{thisseason}.y1,RAD.dataplot{thisseason}.y2,...
            colOBS(thisseason,:),linewidthObs);

        if 1;%thisstation == length(Name)
            ylim([0,0.3]);
            legend(gca,[h(thisseason),hobs(thisseason)],...
                [strcat(Name(thisstation),'-UKCP'),...
                strcat(Name(thisstation),'-Radar')],...
                'Location','NorthEast','fontsize',10);
            legend boxoff
            % legend(gca,[h,hobs],[strcat(Name,'-UKCP'),strcat(Name,'-Radar')]);
            grid minor
            text(500,0.9,getSeasonName(thisseason),'fontsize',12);
        end
    end
end
XYWH = [50,-50,500,500];
set(gcf,'units','points','position',XYWH);

format_xylabel(ha,4,4)

%% Figure 3: Plot two cdfs for given aggregation area from rawIMF
ha = tight_subplot(4,4,[.01 .01],[.1 .01],[.1 .01]);
[h,hobs] = deal([]);
for thisstation = 1:length(Name)
    CPM = load(sprintf('%srawIMF_Area_%s_4Season.mat',dataSP,Name{thisstation}),'rawIMF');
    RAD = load(sprintf('%sObs_rawIMF_Area_%s_4Season.mat',dataSP,Name{thisstation}),'rawIMF');
    for thisseason = 1:4
        thisplot = (thisstation-1)*4+thisseason;%thisseason;
        axes(ha(thisplot));
        aggNo = 1;%24;
        aggArea = CPM.rawIMF{thisseason}(1).areaAgg(aggNo);
        [~,h(thisseason)] = plot_IMFCDF(CPM.rawIMF{thisseason},aggNo,...
            colm(thisstation,:),linewidth);
        
        [~,hobs(thisseason)] = plot_IMFCDF(RAD.rawIMF{thisseason},aggNo,...
            colOBS(thisstation,:),linewidthObs);

        if 1%thisstation == length(Name)
            YL = ylim;
            XL = xlim;
            legend(gca,[h(thisseason),hobs(thisseason)],...
                [strcat(Name(thisstation),'-UKCP'),...
                strcat(Name(thisstation),'-Radar')],...
                'Location','SouthWest','fontsize',10);
            grid minor
            text(XL(1)+0.2,YL(2)*0.50,{sprintf('%4.2f Km^2',aggArea);...
                getSeasonName(thisseason)},'fontsize',10);
        end
    end
end
XYWH = [50,-50,500,500];
set(gcf,'units','points','position',XYWH);
format_xylabel(ha,4,4)



%% AUXILLARY


% percentile 80-99.99 of one location
function [h] = plotPrcRain(E,N,RainEns,scaleF,region,colorV);
%
% Kendon 2002 Fig. 3
%
h = [];

linewidth = 1+1*(numel(RainEns)==1);

colv = colorV;
x = [70,75,77.5,80,85,90,92.5,95,97,99,99.25,99.5,99.7,99.9,99.92,99.95,99.97,99.99];
for enNo = 1:numel(RainEns)
    
    y = zeros(1,numel(x));
    for loci = 1:size(RainEns{enNo},1)
        for locj = 1:size(RainEns{enNo},2)
            rain = double(squeeze(RainEns{enNo}(loci,locj,:)))/scaleF;
            rain = rain(rain>=0.1);
            y = y+prctile(rain,x);
        end
    end
    y = y/50/50;
    
    h = plot(1:numel(x),y,'color',colv(enNo,:),'linewidth',linewidth);
    hold on;
    xlim([2,numel(x)])
    set(gca,'XTick',2:2:numel(x))
    set(gca,'XTickLabel',num2str(x(2:2:end)','%.2f'));
end

end

% radially
function ggg(region,mons,PeriodName,colm,colmOBS,Period)

[y1,y2] = deal([]);

REGIONS = REGIONS_info();
dataSP = 'H:\DATA_CLIMATE\UKCP18\';

linewidth = 1;
linewidthObs = 2;

Pf = [];

[E,N,RainEnsembles,scaleF,region] = loadRainEnsembles(region,mons,Period);
h = plotthis(RainEnsembles,colm,linewidth);

[E,N,RainNimrod,scaleF,region] = loadRainNimrod(region,mons);
hobs = plotthis({RainNimrod},colmOBS,linewidthObs);

legend([h,hobs],{'UKCP18','Radar'},...
    'Location','SouthWest');
ax.YLim = [10^-6,2];
ax.XLim = [2,230];
XLIM = xlim;YLIM = ylim;
text(XLIM(1)*30,YLIM(2)/2,sprintf('%s-%s',region.Name,PeriodName));
xtickformat(gca, '%d');
xticks([1,10,100])
xticklabels({'1','10','100'})


    function h = plotthis(RainEnsembles,colm,linewidth)
        for enNo = 1:length(RainEnsembles)
            colMedian = colm(enNo,:);
            itag = 1;
            R3d = double(RainEnsembles{enNo})/scaleF;
            y1{enNo} = squeeze(R3d(:,:,1));
            y2{enNo} = squeeze(R3d(25,25,:));
            
            %% spatial power spectra
            for iii = 1:size(R3d,3)
                RAPs_tem = squeeze(R3d(:,:,iii));
                RAPs_tem(RAPs_tem<0.1)=0;
                if nanmean(RAPs_tem(:))>0.03 % nanmax(RAPs_tem(:))>20
                    RAPs_tem = RAPs_tem;
                    pl = 0;
                    [f1,Pf{itag},xCell,yCell] = raPsd2d(RAPs_tem,1,pl);
                    itag = itag+1;
                end
            end
            medianPf = nanmean(cell2mat(Pf'),1);
            h = plot_RAPS_space(f1,medianPf,xCell,yCell,colMedian,linewidth);
        end
    end
end


function h = plot_RAPS_space(f1,Pf,xCell,yCell,col,linewidth)

if linewidth >= 2
    linetype = '--';
else
    linetype = '-';
end
h = loglog(2.2*max(f1(:))./f1,Pf,linetype,'LineWidth',linewidth,'color',col);

set(gcf,'color','white')
set(gca,'FontWeight','bold',...
    'XGrid','on','XDir','reverse');

xlabel('Wavelength (km)','FontWeight','Bold');
ylabel('Power','FontWeight','Bold');

% title('Radially averaged power
% spectrum','FontSize',fontSize,'FontWeight','Bold');
hold on;

end

function [h1,h2] = plot_WetArea(Thres,y1,y2,colv,lineType,lineWidth)

setFigureProperty('Subplot2')
[h1,h2] = deal(NaN);

%     subplot(1,2,1);
h1 = plot(repmat(Thres,size(y1,1),1)',y1',lineType,'color',colv,...
    'Markerfacecolor','w','Markersize',3,'Linewidth',lineWidth);
h1 = h1(1);
configPlot_war()
ylabel('WAR','Interpreter','Latex')
hold on;

ytickformat(gca, '%g%%');
%     subplot(1,2,2);
%     h2 = plot(repmat(Thres,12,1),y2,'.-','color',colv);
%     h2 = h2(1);
%     configPlot_war()
%     ylabel('$std(\bar{A}(R<Thre))$','Interpreter','Latex')
%     hold on;


end


function [h1,h2] = plot_IMF(areaAgg,y1,y2,colv,linewidth)

setFigureProperty('Subplot2')
[h1,h2] = deal([]);
    
%     h2 = plot(repmat(areaAgg,size(y2,1),1)',y1','-','color',colv,'Linewidth',linewidth);
%     h2 = h2(1);
%     configPlot_imf()
%     ylabel('Mean of IMF')
%     hold on;
    
%%%%%%%%%%
    h2 = plot(repmat(areaAgg,size(y2,1),1)',y2','-','color',colv,'Linewidth',linewidth);
    h2 = h2(1);
    configPlot_imf()
    ylabel('Mean of IMF(IMF>0)')
    hold on;

    
end


function [h1,h2] = plot_IMFCDF(rawIMF,aggNo,colv,linewidth)

h1 = 0;

for enNo = 1:length(rawIMF)
    [f,x,flo,fup] = ecdf(rawIMF(enNo).AreaIMFs{aggNo});
    h2 = plot(x,1-f,'color',colv,'linewidth',linewidth);
    hold on;
end

configPlot_imfcdf()

end


function configPlot_imfcdf()
set(gca,'YScale','log')
XL = xlim;
xlim([2,XL(2)])
grid minor
YL = ylim;
ylim([10^-5,YL(2)])
% legend([h,hobs],{'R1-UKCP','R1-Radar'});
xlabel('Areal Mean Rainrate mm/h')
ylabel('Exceedence probability');
end


function configPlot_war()

% set(gca, 'YScale', 'log')
xlabel('Threshold [mm/h]');
grid minor

end

function configPlot_imf()

xlabel('Area [km^2]');
ylabel('mean rain rate mm/h');
% set(gca, 'XScale', 'log')
grid minor

end


