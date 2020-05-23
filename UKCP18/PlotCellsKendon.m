clear;clc
close all

addpath(genpath('C:\Users\Yuting Chen\Dropbox (Personal)\MatLab_Func'));
addpath(genpath(cd));


REGIONS = REGIONS_info();

region = REGIONS.London;%Westuk;%SWestuk;%Scotland;%


plotOneRegionCellDensityHist(REGIONS.London)
plotOneRegionCellDensityHist(REGIONS.Westuk)
plotOneRegionCellDensityHist(REGIONS.SWestuk)
plotOneRegionCellDensityHist(REGIONS.Scotland)


function plotOneRegionCellDensityHist(region)

dataSP = 'H:\DATA_CLIMATE\UKCP18\';

% figure;
% setFigureProperty('Meeting');
% XYWH = [50,-50,500,360];
% set(gcf,'units','points','position',XYWH);
% ha = tight_subplot(2,2,[.08 .08],[.1 .08],[.1 .08]);


figure;
setFigureProperty('Meeting');
XYWH = [50,-50,500,360];
set(gcf,'units','points','position',XYWH);
ha2 = tight_subplot(2,4,[.18 .06],[.1 .08],[.1 .08]);


figure
ha3 = tight_subplot(2,2,[.08 .08],[.15 .05],[.15 .05]);

colm = pink(30);
linewidth = 1.5;
linewidthObs = 3;

f1 = @(x)sqrt(x);
f2 = @(x)x;%ceil(x);

STAT = struct;
STAT.SIMNUM = [];
STAT.OBSNUM = [];
pind = 0;
for seasonNo = [2,3,4,1]%1:4
    pind = pind+1;
    % legend(h,)
    OBS = load(sprintf('%sCellProp_NIMROD_Season%01d_%s.mat',dataSP,seasonNo,region.Name),...
        'CELLA','THRE');
    CPM = load(sprintf('%sCellProp_UKCP_Season%01d_%s.mat',dataSP,seasonNo,region.Name),...
        'CELLA','THRE');
    Yobs = [];Ysim = [];
    
    for i = 1:length(OBS.THRE)-1% exclude the last one (100m/h)
        
        Yobs{i} = OBS.CELLA{i}{1};
        
        Yobs{i} = f2(f1(Yobs{i}));
    end
    
    for i = 1:length(CPM.THRE)-1% exclude the last one (100m/h)
        Ysim{i} = [];
        for enNo = 1:12
            Ysim{i} = [Ysim{i};CPM.CELLA{i}{enNo}];
        end
        Ysim{i} = f2(f1(Ysim{i}));
    end
    
    CPM.THRE = CPM.THRE(1:end-1);
    OBS.THRE = OBS.THRE(1:end-1);
    %% FIGURE 1: PLOT VIOLIN
    %{
    axes(ha(seasonNo))
    grid minor
    text(6.5,25,sprintf('%s-%s',region.Name,getSeasonName(seasonNo)),'fontsize',12,...
        'HorizontalAlignment','right');
    hansim = distributionPlot(Ysim,'widthDiv',[2 1],'histOri','left',...
        'color',[0.8588,0.4392,0.5765],'histOpt',0,'showMM',4);%,'divFactor',50);%,'globalNorm',3)
    hanobs = distributionPlot(ha(seasonNo),Yobs,'widthDiv',[2 2],'colormap',[0.4667,0.5333,0.6000],...
        'histOri','right',...
        'histOpt',0,'showMM',4,...%,'divFactor',1/1.5,...
        'xNames',num2cell(OBS.THRE),'globalNorm',0);
    EDGES = sqrt(([5,20,50,100,200,500]-0.01));%-0.05;
    set(gca,'YTick',round(EDGES));
    set(gca,'YTickLabel',2.2^2*round((EDGES+0.01).^2));
    % ylim([-1 50])%'colormap',copper,
    ylabel('Cell size (Km^2)');
    legend([hanobs{1}{1},hansim{1}{1}],{'Radar','CPMs'},'Location','NorthWest')
    set(ha(seasonNo), 'YAxisLocation', 'left')
    % cptcmap('GMT_polar', 'mapping', 'scaled');
    box on;
    %}
    
    %% FIGURE 2: PLOT KENDON 2012 Figures;
    CELLSIZE = [1,5,10,20,50,100,200,500,2500,inf];
    EDGES = sqrt((CELLSIZE-0.01));%-0.05;
    axes(ha2(seasonNo*2-1))
    [X,Xtrue,Y,GMAPsim] = plotKendon2012(Ysim,CPM.THRE,EDGES,CELLSIZE);
    title(sprintf('%s-CPM 2.2Km (1980-2000)',getSeasonName(seasonNo)));
    
    axes(ha2(seasonNo*2))
    [X,Xtrue,Y,GMAPobs] = plotKendon2012(Yobs,CPM.THRE,EDGES,CELLSIZE);
    title(sprintf('%s-UK Weather Radar (2007-2018)',getSeasonName(seasonNo)));
    
    axes(ha3(pind))
    plotKendonBias(X,Xtrue,Y,GMAPsim-GMAPobs);
    title(sprintf('(%s)CPM-Radar',getSeasonName(seasonNo)));
    
    STAT.SIMNUM{seasonNo} = cell2mat(cellfun(@(x)length(x)/12/20, Ysim, 'UniformOutput', false));
    STAT.OBSNUM{seasonNo} = cell2mat(cellfun(@(x)length(x)/11, Yobs, 'UniformOutput', false));
    
end


format_xylabel(ha3,2,2)
filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\CellProperty';
filename = [filePath,filesep,'cellPD_',region.Name];
savePlot(filename,'XYWH',[150,0,600,500],'needreply','N');

%%
figure;
ha4 = tight_subplot(2,2,[.15 .05],[.15 .1],[.15 .05]);
pind = 0;
for seasonNo = [2,3,4,1]
    pind = pind+1;
    axes(ha4(pind));
    plot(CPM.THRE,STAT.SIMNUM{seasonNo},'ro-','linewidth',2);
    hold on;
    plot(CPM.THRE,STAT.OBSNUM{seasonNo},'ko-','linewidth',2);
    title(sprintf('%s',getSeasonName(seasonNo)))
    ylabel('cell nums/year');
    xlabel('threshold mm/h')
    legend('SIM','OBS')
    % set(gca,'YScale','log');
    set(gca,'linewidth',2);
    grid minor
    xlim([0,20])
    ylim([0,4000])
end
format_xylabel(ha4,2,2)

filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\CellProperty';
filename = [filePath,filesep,'cellNumbers_',region.Name];
savePlot(filename,'XYWH',[150,0,600,500],'needreply','N');

close all
end



%%
function [X,Xtrue,Y,GMAP] = plotKendon2012(Z,THRE,EDGES,CELLSIZE)

GMAP = NaN(length(THRE),numel(EDGES)-1);

for i = 1:length(THRE)
    Z{i}(Z{i}<=0)=[];
    [H,~] =  histcounts(Z{i},EDGES);
    GMAP(i,:) = H./sum(H);
    %[GMAP(i,:),EDGES] =  histcounts(Ysim{i});
end

GMAP(GMAP==0)=NaN;



X = 1:length(EDGES(1:end-1));
Xtrue = CELLSIZE(1:end-1);
Y = THRE;
%{
pcolor(X,Y,GMAP);shading flat

[XX,YY] = meshgrid(X,Y);

GMAP(isnan(GMAP))=0;

fxshow = @(x)reshape(x(2:end,1:end-1),[],1);
text(fxshow(XX),fxshow(YY)*0.9,sprintfc('%.0f%%',round(100*fxshow(GMAP))),...
    'fontsize',9,'fontweight','bold','HorizontalAlignment','left',...
    'VerticalAlignment','top')

cptcmap('GMT_haxby', 'mapping', 'direct','ncol',20,'flip',true);
set(gca,'clim',[0,0.5]);
box on
xlim([1-0.01 length(EDGES(1:end-1))+0.01]);
ylim([1-0.01,20+0.01])
set(gca,'XTick',1:length(EDGES(1:end-1)))
set(gca,'XTickLabel',CELLSIZE(1:end-1));
xlabel('cell size ((2.2km)^2 grid pts)')
ylabel('threshold (mm/h)')
colorbar
%}

end

function plotKendonBias(X,Xtrue,Y,Bias)

Y0 = 1:length(Y);
pcolor(X,Y0,Bias);% shading flat

[XX,YY] = meshgrid(X,Y);

% GMAP(isnan(GMAP))=0;

% fxshow = @(x)reshape(x(2:end,1:end-1),[],1);
% text(fxshow(XX),fxshow(YY)*0.9,sprintfc('%.0f%%',round(100*fxshow(GMAP))),...
%     'fontsize',9,'fontweight','bold','HorizontalAlignment','left',...
%     'VerticalAlignment','top')

cptcmap('diff_4_bias', 'mapping', 'scaled','ncol',9);
set(gca,'clim',[-0.3,0.3]);
box on
% xlim([1-0.01 length(EDGES(1:end-1))+0.01]);
% ylim([1-0.01,20+0.01])
set(gca,'YTick',Y0)
set(gca,'YTickLabel',Y);
set(gca,'XTick',X)
set(gca,'XTickLabel',Xtrue);
xlabel('cell size ((2.2km)^2 grid pts)')
ylabel('threshold (mm/h)')
colorbar
ax = gca;
ax.XTickLabelRotation = 90;
end






