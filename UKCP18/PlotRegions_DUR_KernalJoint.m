close all
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



%% Figure 1: Plot the saved dataplot

[h,hobs] = deal([]);

logsf = 1;%5;


for thisstation = 1:length(Name)
    figure;
    setFigureProperty('Meeting');
    XYWH = [50,-50,500,500];
    set(gcf,'units','points','position',XYWH);
    ha = tight_subplot(2,4,[.08 .05],[.15 .1],[.08 .12]);
    cmap = pink(length(Name)*2);
    
    figure;
    setFigureProperty('Meeting');
    ha_bias = tight_subplot(2,2,[.10 .1],[.15 .1],[.1 .14]);
    tsp=1;
    
    for thisseason = [2,3,4,1]%[2 3 4 1]
        thisplot = thisseason;% different location in same figure;%(thisstation-1)*4+tsp;%thisseason;
        
        CPM = load(sprintf('%sDUR_Area_%s_1980-2000_4Season.mat',dataSP,Name{thisstation}),'dataplot');
        RAD = load(sprintf('%sObs_DUR_Area_%s_4Season.mat',dataSP,Name{thisstation}),'dataplot');
        
        ThreDur = [1:24,30,36,42,48,60,72];%[0,1,3,6,12,24,48,100,1000];%[1:2:5,6:10];
        ThreImf = [0.1 0.2 0.5 1 2 5 10 20];%[0.03,0.1,0.2,0.4,0.6,1,2,3,5];%[0:0.1:1,1.5,2,2.5,3];

        
        % do it for CPM
        axes(ha(tsp));
        Dur1 = arrayfun(@(x)cell2mat(cellfun(@(y)y.Dur,x,'UniformOutput', false)),...
            CPM.dataplot{thisseason},'UniformOutput', false);
        Imf1 = arrayfun(@(x)cell2mat(cellfun(@(y)y.pIMF,x,'UniformOutput', false)),...
            CPM.dataplot{thisseason},'UniformOutput', false);
        
        [X,Y,GMAP_cpm] = plotJointDist([cell2mat(Imf1'),cell2mat(Dur1')],ThreImf,ThreDur);
        % plotBoxWAR(M,THRE,EDGES);
        formatIt('CPM2.2',thisseason)
        
        % do it for radar
        axes(ha(tsp+4));
        [X,Y,GMAP_rad] = plotJointDist([RAD.dataplot{thisseason}.pIMF,...
            RAD.dataplot{thisseason}.Dur],ThreImf,ThreDur);
        formatIt('Radar',thisseason)
        
        % do it for bias
        axes(ha_bias(tsp));
        contourf(log2(1+X),log2(Y),GMAP_cpm-GMAP_rad,[-0.06:0.01:0.06]);shading flat
        formatIt('Radar',thisseason)
        cptcmap('diff_4_bias', 'mapping', 'scaled','ncol',19);
        set(gca,'clim',[-0.04,0.04]);
        box on
        % ylim([1-0.01,20+0.01])
        set(gca,'YTick',log2([0.1 0.2 0.5 1 2 5 10 20]));
        set(gca,'YTickLabel',[0.1 0.2 0.5 1 2 5 10 20]);
        set(gca,'XTick',log2(1+[1,5,10,15,20,24,36,48,96]));
        set(gca,'XTickLabel',[1,5,10,15,20,24,36,48,96]);

        xlabel('Duration (hour)')
        ylabel('Peak Intensity (mm/h)')
        formatIt('CPM-Radar',thisseason)
        
        tsp=tsp+1;
    end
    
    
    XYWH = [50,-50,500,500];
    set(gcf,'units','points','position',XYWH);
    
    axes(ha(end))
    format_xylabel(ha,2,4)
    c = colorbar('location','Manual', 'position', [0.94 0.1 0.015 0.8]);
    c.Ruler.TickLabelFormat='%0.4f';%'%g%%';
    
    filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\TemporalAnalysis\';
    filename = [filePath,filesep,'Kernel_DUR_PI_rad_cpm_hist3_Allvalue_',Name{thisstation}];
    % savePlot(filename,'XYWH',[150,0,1400,700],'needreply','N');
    
    axes(ha_bias(end))
    format_xylabel(ha_bias,2,2)
    c = colorbar('location','Manual', 'position', [0.92 0.12 0.015 0.78]);
    c.Ruler.TickLabelFormat='%g';
    
    filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\TemporalAnalysis\';
    filename = [filePath,filesep,'Kernel_DUR_PI_rad_cpm_hist3_',Name{thisstation}];
    % savePlot(filename,'XYWH',[150,0,700,800],'needreply','N');
    pause(1)
    close all
    
end


function formatIt(productName,thisseason)
set(gca,'linewidth',1)
box on
title(sprintf('%s (%s)',productName,getSeasonName(thisseason)));
% text(0,0.999,getSeasonName(thisseason),'fontsize',12);
ax = gca;
% ax.FontSize = 12;
ax.XTickLabelRotation = 90;
set(gca,'TickDir','out')
xlim([ax.XLim(1)-0.1 ax.XLim(2)+0.1]);
ylim([ax.YLim(1)-0.1,ax.YLim(2)+0.1])
end


function [X,Y,GMAP] = plotJointDist(XY,ThreX,ThreY)

% figure;
x = log2(XY(:,2));y = log2(XY(:,1));
% Estimate a continuous pdf from the discrete data
gridx1 = linspace(log2(1),log2(100),100);%linspace(0,7,100);% linspace(min(x),max(x),50);
gridx2 = linspace(-3,7,100);% linspace(min(y),max(y),100);
[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
hxy = ksdensity([x,y],xi,'BoundaryCorrection','reflection','PlotFcn','surf','Kernel','epanechnikov');
hxy = reshape(hxy,numel(gridx2),numel(gridx1));
surf(gridx1,gridx2,hxy);
contourf(gridx1,gridx2,hxy,[0:0.02:0.12]);shading interp;
hxy(hxy==0) = NaN;
shading flat
set(gca,'XLim',[min(x) max(x)])
set(gca,'YLim',[min(y) max(y)])
% mesh(xxi,yyi,pdfxy)
% contourf(xxi,yyi,pdfxy);shading interp;
cptcmap('flood_blue','mapping','scaled','ncol',100);
caxis([0,0.12]);
xlabel('Duration')
ylabel('Intensity')
ax = gca;
ax.XTickLabel = 2.^ax.XTick;
ax.YTickLabel = 2.^ax.YTick;
axis on
box on

X = 2.^gridx1;
Y = 2.^gridx2;
GMAP = hxy;


box on
xlabel('Duration (hour)')
ylabel('Peak Intensity (mm/h)')

end



