                                         
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
    ha = tight_subplot(2,4,[.08 .05],[.15 .1],[.1 .1]);
    cmap = pink(length(Name)*2);
    
    figure;
    setFigureProperty('Meeting');
    ha_bias = tight_subplot(2,2,[.10 .1],[.15 .1],[.1 .14]);
    tsp=1;
    
    for thisseason = [2,3,4,1]%[2 3 4 1]
        thisplot = thisseason;% different location in same figure;%(thisstation-1)*4+tsp;%thisseason;
        
        CPM = load(sprintf('%sWAR_Area_%s_1980-2000_4Season.mat',dataSP,Name{thisstation}),'dataplot');
        RAD = load(sprintf('%sObs_WAR_Area_%s_4Season.mat',dataSP,Name{thisstation}),'dataplot');
        
        thre_ind = [1:2:5,6:10];
        
        thre = CPM.dataplot{1, 1}(1).Thres;
        THRE = thre(thre_ind);%[1,5,10,20,50,100,200,500,2500,inf];
        EDGES = [0,5e-3,1e-2,5e-2,1e-1,0.5,0.7,1,Inf];%THRE;%sqrt((CELLSIZE-0.01));%-0.05;
        % [0,5e-3,1e-2:0.01:5e-2,1e-1:0.1:0.5,0.7:0.1:1,Inf]% 
        % do it for CPM
        axes(ha(tsp));
        M = arrayfun(@(x)x.y2(thre_ind),CPM.dataplot{thisseason},...
            'UniformOutput', false);
        M = mat2cell(cell2mat(transpose(cat(1,M{:}))),ones(1,length(THRE)));
        [X,Y,GMAP_cpm] = plotBoxWAR(M,THRE,EDGES);
        formatIt('CPM2.2',thisseason)
        
        % do it for radar
        axes(ha(tsp+4));
        [X,Y,GMAP_rad] = plotBoxWAR(RAD.dataplot{thisseason}(1).y2(thre_ind),THRE,EDGES);
        formatIt('Radar',thisseason)

        % do it for bias
        axes(ha_bias(tsp));
        pcolor(X,Y,100*(GMAP_cpm-GMAP_rad)./GMAP_rad);% shading flat
        formatIt('Radar',thisseason)
        cptcmap('diff_4_bias', 'mapping', 'scaled','ncol',19);
        set(gca,'clim',[-200,200]);
        box on
        xlim([1-0.01 length(EDGES(1:end-1))+0.01]);
        % ylim([1-0.01,20+0.01])
        set(gca,'XTick',1:length(EDGES(1:end-1)))
        set(gca,'XTickLabel',EDGES(1:end-1));
        xlabel('Wet Area Ratio')
        ylabel('threshold (mm/h)')
        formatIt('CPM-Radar',thisseason)

        tsp=tsp+1;
    end
    

XYWH = [50,-50,500,500];
set(gcf,'units','points','position',XYWH);

axes(ha(end))
format_xylabel(ha,2,4)
c = colorbar('location','Manual', 'position', [0.92 0.12 0.015 0.8]);
c.Ruler.TickLabelFormat='%g%%';

filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\SpatialAnalysis\';
filename = [filePath,filesep,'WAR_rad_cpm_histCountAllvalue_',Name{thisstation}];
savePlot(filename,'XYWH',[150,0,1400,700],'needreply','N');

axes(ha_bias(end))
format_xylabel(ha_bias,2,2)
c = colorbar('location','Manual', 'position', [0.9 0.12 0.015 0.78]);
c.Ruler.TickLabelFormat='%g%%';

filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\';
filename = [filePath,filesep,'WAR_rad_cpm_histCount_',Name{thisstation}];
savePlot(filename,'XYWH',[150,0,700,800],'needreply','N');

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
% set(gca,'linewidth',2)
set(gca,'TickDir','out')
% xlim([ax.XLim(1)-0.1 ax.XLim(2)+0.1]);
% ylim([ax.YLim(1)-0.1,ax.YLim(2)+0.1])
end


function [X,Y,GMAP] = plotBoxWAR(Z,THRE,EDGES)

GMAP = NaN(length(THRE),numel(EDGES)-1);

for i = 1:length(THRE)
    Z{i}(Z{i}<=0)=[];
    [H,~] =  histcounts(Z{i},EDGES);
    GMAP(i,:) = H./sum(H);
    %[GMAP(i,:),EDGES] =  histcounts(Ysim{i});
end

GMAP(GMAP==0)=NaN;

X = 1:length(EDGES(1:end-1));
Y = THRE;
pcolor(X,Y,GMAP);% shading flat

% [XX,YY] = meshgrid(X,(Y(1:end-1)+Y(2:end))/2);
% 
% GMAP(isnan(GMAP))=0;
% 
% fxshow = @(x)reshape(x(2:end,1:end-1),[],1);
% text(fxshow(XX),fxshow(YY),sprintfc('%.0f%%',round(100*fxshow(GMAP))),...
%     'fontsize',9,'fontweight','bold','HorizontalAlignment','left',...
%     'VerticalAlignment','top')

cptcmap('GMT_haxby', 'mapping', 'scaled','ncol',10,'flip',true);
set(gca,'clim',[0,1]);
box on
xlim([1-0.01 length(EDGES(1:end-1))+0.01]);
% ylim([1-0.01,20+0.01])
set(gca,'XTick',1:length(EDGES(1:end-1)))
set(gca,'XTickLabel',EDGES(1:end-1));
xlabel('Wet Area Ratio (110km)')
ylabel('threshold (mm/h)')

end



