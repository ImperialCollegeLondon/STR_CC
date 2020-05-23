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




for thisstation = 1:length(Name)
    subplot(1,4,thisstation)
    
    CPM = load(sprintf('%sDUR_Area_%s_1980-2000_4Season.mat',dataSP,Name{thisstation}),'dataplot');
    RAD = load(sprintf('%sObs_DUR_Area_%s_4Season.mat',dataSP,Name{thisstation}),'dataplot');
    
    for ensNo = 1:12
        A(ensNo,:) = cellfun(@(c)length(c{ensNo}.stormsi)/20/3,CPM.dataplot);
        
        hold on;
    end
    plot(1:4,nanmedian(A,1),'r-');
    
    hsimran = fill([1:4,4:-1:1],[nanmin(A,[],1),flip(nanmax(A,[],1))],'r',...
        'LineStyle','none');hold on;
    alpha(0.3)
    
    A = cellfun(@(c)length(c.stormsi)/12/3,RAD.dataplot);
    plot(1:4,A,'k-')
    hold on
    
    title(Name{thisstation})
    ylabel('storm/month');
    box on
    set(gca,'linewidth',1)
    ax = gca;
    ax.XTick = 1:4;
    ax.XTickLabel = arrayfun(@(s)getSeasonName(s),1:4,'UniformOutput',false);
    ax.YLim=[20,30];%ax.YLim(2)];
end