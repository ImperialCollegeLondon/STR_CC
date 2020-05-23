clear;clc

addpath(genpath('C:\Users\Yuting Chen\Dropbox (Personal)\MatLab_Func'));
addpath(genpath(cd));

REGIONS = REGIONS_info();
region = REGIONS.UK;

%% DENSITY PLOT plot+save <cell.centroid> for whole uk

figure;
ha_bias = tight_subplot(2,2,[.05 .05],[.02 .02],[.02 .10]);


upperl = [];
[CELLALL_rad,CELLALL_cpm] = deal(table);
seasons = 2;
for season = seasons
    
    % LOAD
    
    source = 'RAD';
    [E_rad,N_rad,CELLALL_rad0,THRE] = getAllCellsUK(season,region,source,1);
    
    source = 'CPM';
    [E_cpm,N_cpm,CELLALL_cpm0,THRE] = getAllCellsUK(season,region,source,1);
    
    CELLALL_rad = [CELLALL_rad;CELLALL_rad0];
    CELLALL_cpm = [CELLALL_cpm;CELLALL_cpm0];
    
end
% PLOT & SAVE
%
figure;
setFigureProperty('Paper');
ha = tight_subplot(1,3,[.05 .05],[.05 .05],[.05 .10]);
set(gcf,'units','centimeters','position',[5,0,24,11])
plotComp = true;

axes(ha(1))
source = 'RAD';
[E,N,Count_1,Count_obs] = plotCells(E_rad,N_rad,CELLALL_rad,seasons,region,source,'SameArea',plotComp);%'TRIM');
axis equal;xlim([-50,750]);ylim([0,1300])
caxis([0,2])
axes(ha(2))
source = 'CPM';
[E,N,Count_2,Count_sim] = plotCells(E_cpm,N_cpm,CELLALL_cpm,seasons,region,source,'SameArea',plotComp);%'TRIM');
axis equal;xlim([-50,750]);ylim([0,1300])
c = colorbar(ha(2),'location','Manual', 'position', [0.55,0.57,0.02,0.22]);
c.Ticks = [0,2];
c.TickLabels(end) = strcat(c.TickLabels(end),'/year');


sig = testSig(Count_sim,Count_obs);



axes(ha(3))
f = @(mat)imresize(mat,0.2,'box');%mat;%imresize(mat,0.2,'box');
BIA = (f(Count_2)-f(Count_1))./f(Count_1)*100;
% BIA(sig==0) = 0;
pcolor(f(E),f(N),BIA);
hold on;
shading flat
if isempty(upperl)
    upperl = 150;%prctile(Count_2(Count_2>0),90);
end
caxis(gca,[-upperl,upperl])
cptcmap('diff_darkBlue_darkRed','mapping', 'scaled');%,'ncol',upperl*10+1);
% alpha(0.8);
xlabel('Easting/Km');
UKMap = load('H:\CODE_MATLAB\SpatialTemporalDATA\shapeFileFolder\UKBorderGrid.mat');
plot(UKMap.borderE/1000,UKMap.borderN/1000,'k-');
ylabel('Northing/Km')
XYWH = [150,0,420,490];
set(gcf,'units','points','position',XYWH);
axis equal
text(-50,1300,sprintf('%%Difference (CPM-Radar)'),...
    'fontweight','bold','fontsize',14,...
    'verticalalignment','top',...
    'horizontalalignment','left');

% colorbar;
xlim([-50,750]);ylim([0,1300])
axis off
c = colorbar(ha(3),'location','Manual', 'position', [0.85,0.57,0.02,0.22]);
c.Ticks = [-100,0,100];
c.TickLabels = strcat(c.TickLabels,'%');

if plotComp
    dataPath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP';
    if numel(seasons) == 1
        filename = [dataPath,filesep,'CellDensity_',getSeasonName(seasons),'_Rad_Cpm',sprintf('_%02dmm',THRE)];
    else
        filename = [dataPath,filesep,'CellDensity_Rad_Cpm',sprintf('_%02dmm',THRE)];
    end
    savePlot(filename,'units','centimeters','XYWH',[5,0,24,11],'needreply','Y','onlyPng',true);
    % close ha
end

close all



%%
%     axes(ha_bias(season))
%     f = @(mat)imresize(mat,0.2,'box');
%     pcolor(f(E),f(N),f(Count_2)-f(Count_1));
%     hold on;
%     shading flat
%     if isempty(upperl)
%     upperl = prctile(Count_2(Count_2>0),90);
%     end
%     caxis(gca,[-upperl,upperl])
%     cptcmap('diff_darkBlue_darkRed','mapping', 'scaled','ncol',upperl*10);
%     % alpha(0.8);
%     xlabel('Easting/Km');
%     UKMap = load('H:\CODE_MATLAB\SpatialTemporalDATA\shapeFileFolder\UKBorderGrid.mat');
%     plot(UKMap.borderE/1000,UKMap.borderN/1000,'k-');
%     ylabel('Northing/Km')
%     XYWH = [150,0,420,490];
%     set(gcf,'units','points','position',XYWH);
%
%     text(-50,1300,sprintf('Difference CPM-Radar (%s)',getSeasonName(season)),...
%         'fontweight','bold','fontsize',14,...
%         'verticalalignment','top',...
%         'horizontalalignment','left');
%
%     % colorbar;
%     xlim([-50,750]);
%     ylim([0,1300])
%     axis off



% colorbar('location','Manual', 'position', [0.93 0.1 0.02 0.81]);

%
% dataPath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP';
% filename = [dataPath,filesep,'CellDensity_4Season_bias_Rad_Cpm',sprintf('_%02dmm',THRE)];
% savePlot(filename,'XYWH',[150,0,600,800]);%,'needreply','N');


%% HIST plot+save <cell.area> for each region ALL <region> in <REGIONS>;

for season = 2%
    source = 'RAD';
    [E_rad,N_rad,CELLALL_rad,THRE] = getAllCellsUK(season,REGIONS.UK,source,2);
    
    source = 'CPM';
    [E_cpm,N_cpm,CELLALL_cpm,THRE] = getAllCellsUK(season,REGIONS.UK,source,2);
    
    region = REGIONS.SCO;% SWestuk;%London;
    try
        REGIONS = rmfield(REGIONS,{'smallUK','Birm','London','EAng','Scotland','Westuk','SWestuk'});
    catch
    end
    shortNames = structfun(@(region)plotAreaHist(region,CELLALL_rad,CELLALL_cpm,season,THRE), ...
        REGIONS, 'UniformOutput', false);
    
    shortNames = structfun(@(region)plotEccencHist(region,CELLALL_rad,CELLALL_cpm,season,THRE), ...
        REGIONS, 'UniformOutput', false);
    close all
end

%%
load('H:\CODE_MATLAB\SpatialTemporalDATA\shapeFileFolder\UKBorderGrid.mat')

load(['K:\UK_shape\DTM50.mat'],'DTM50','Eno','Nno')

DTM_2km = imresize(DTM50,0.025);
clear DTM50
Eno = Eno(1):2000:Eno(end);
Nno = Nno(1):-2000:Nno(end);
[EE,NN] = meshgrid(Eno,Nno);
in = inpolygon(EE/1000,NN/1000,borderE/1000,borderN/1000);
DTM_2km(~in) = NaN;

%%

figure;

setFigureProperty('Paper')

ax1 = axes;
s = pcolor(Eno/1000,Nno/1000,DTM_2km);
shading flat
hold on;

view(2)
ax2 = axes;
Count = histcounts2(CELLALL.Centroid(:,2),CELLALL.Centroid(:,1),size(LAT));
Count = Count./nansum(Count(:));
h = pcolor(ax2,E,N,Count);
alpha(h,0.25)
hold on;
shading flat


%%Link them together
linkaxes([ax1,ax2])
%%Hide the top axes
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
%%Give each one its own colormap
caxis(ax2,[prctile(Count(:),20),prctile(Count(:),80)])
cptcmap('GMT_no_green',ax2, 'mapping', 'scaled')
cptcmap('GMT_gray',ax1, 'mapping', 'scaled','flip',true);


hold on;
plot(borderE/1000,borderN/1000,'-','color',ones(1,3)*0.5,'linewidth',1);


%%
%%Then add colorbars and get everything lined up
set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb1 = colorbar(ax1,'Position',[.08 .11 .0275 .815]);
cb2 = colorbar(ax2,'Position',[.88 .11 .0275 .815]);

ylim([-200,1200])
xlim([-200,800])
text(-180,1150,sprintf('Rain Cell Density Map - %s (2.2km)',getSeasonName(season)),...
    'fontweight','bold','fontsize',12);

%%
filename = [cd,filesep,'CellDensity_MAM_CPM'];
savePlot(filename,'XYWH',[150,0,600,700]);

%% AUXILLARY FUNC FOR PLOT
function sig = testSig(Count_sim,Count_obs)
sig = ttest(permute(Count_sim-Count_obs,[3,1,2]),0,'alpha',0.05);
sig = squeeze(sig);
end
function ha = plotEccencHist(region,CELLALL_rad,CELLALL_cpm,season,thre)

borderE = [region.minE, region.minE, region.minE+region.dx*(region.dimE-1),region.minE+region.dx*(region.dimE-1)];
borderN = [region.minN, region.minN+region.dx*(region.dimN-1), region.minN+region.dx*(region.dimN-1),region.minN];

[in] = inpolygon(CELLALL_rad.Centroid(:,1)*2.2-280,CELLALL_rad.Centroid(:,2)*2.2-226,borderE,borderN);
rmInd = CELLALL_rad.Area==1 | CELLALL_rad.Area>1000000 | (~in);
Y_rad = CELLALL_rad.Eccentricity(~rmInd); %Area(~rmInd); %

[in] = inpolygon(CELLALL_cpm.Centroid(:,1)*2.2-312,CELLALL_cpm.Centroid(:,2)*2.2-229,borderE,borderN);
rmInd = CELLALL_cpm.Area==1 | CELLALL_cpm.Area>1000000 | (~in);
Y_cpm = CELLALL_cpm.Eccentricity(~rmInd); %Area(~rmInd); %

figure;
ha = tight_subplot(1,2,[.05 .05],[.25 .05],[.15 .05]);

axes(ha(1))
histogram(Y_rad,'BinWidth',0.05,'Normalization','probability');
text(0.9,0.5,sprintf('RAD-%s(%s)',region.Name,getSeasonName(season)),'horizontalalignment','right',...
    'verticalalignment','top')
formatHistEccenc()

axes(ha(2))
h = histogram(Y_cpm,'BinWidth',0.05,'Normalization','probability');
text(0.9,0.5,sprintf('CPM-%s(%s)',region.Name,getSeasonName(season)),'horizontalalignment','right',...
    'verticalalignment','top')
formatHistEccenc()

format_xylabel(ha,1,2)


dataPath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\CellProperty';
filename = [dataPath,filesep,sprintf('Cell_Eccenc_Hist_Rad_Cpm_%s_%s_thre%03dmm',region.Name,getSeasonName(season),thre)];
savePlot(filename,'XYWH',[150,100,400,200],'needreply','N','onlyPng',true);

close all
    function formatHistEccenc()
        xlim([0,1]); ylim([0,0.5])
        ylabel('pdf')
        xlabel('Eccentricity');
        grid minor
    end
end



function ha = plotAreaHist(region,CELLALL_rad,CELLALL_cpm,season,thre)

borderE = [region.minE, region.minE, region.minE+region.dx*(region.dimE-1),region.minE+region.dx*(region.dimE-1)];
borderN = [region.minN, region.minN+region.dx*(region.dimN-1), region.minN+region.dx*(region.dimN-1),region.minN];

[in] = inpolygon(CELLALL_rad.Centroid(:,1)*2.2-280,CELLALL_rad.Centroid(:,2)*2.2-226,borderE,borderN);
rmInd = CELLALL_rad.Area==1 | CELLALL_rad.Area>1000000 | (~in);
Y_rad = CELLALL_rad.Area(~rmInd); % Eccentricity(~rmInd); %

[in] = inpolygon(CELLALL_cpm.Centroid(:,1)*2.2-312,CELLALL_cpm.Centroid(:,2)*2.2-229,borderE,borderN);
rmInd = CELLALL_cpm.Area==1 | CELLALL_cpm.Area>1000000 | (~in);
Y_cpm = CELLALL_cpm.Area(~rmInd); % Eccentricity(~rmInd); %

figure;
ha = tight_subplot(1,2,[.05 .05],[.25 .05],[.15 .05]);

axes(ha(1))
histogram(Y_rad,'BinWidth',5,'Normalization','probability');
text(190,0.5,sprintf('RAD-%s(%s)',region.Name,getSeasonName(season)),'horizontalalignment','right',...
    'verticalalignment','top')
formatHistArea()

axes(ha(2))
h = histogram(Y_cpm,'BinWidth',5,'Normalization','probability');
text(190,0.5,sprintf('CPM-%s(%s)',region.Name,getSeasonName(season)),'horizontalalignment','right',...
    'verticalalignment','top')
formatHistArea()

format_xylabel(ha,1,2)


dataPath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\CellProperty';
filename = [dataPath,filesep,sprintf('Cell_Area_Hist_Rad_Cpm_%s_%s_thre%03dmm',region.Name,getSeasonName(season),thre)];
savePlot(filename,'XYWH',[150,100,400,200],'needreply','N','onlyPng',true);

close all

    function formatHistArea()
        xlim([0,199]); ylim([0,0.5])
        ylabel('pdf')
        xlabel('area(2.2km grid pts)');
        grid minor
    end
end



function [E,N,Count,CountAll] = plotCells(E,N,CELLALL,seasons,region,source,option,pl)

arguments
    E (:,:) double
    N (:,:) double
    CELLALL (:,:) table
    seasons (1,:) double
    region (1,1) struct
    source (1,:) char {mustBeMember(source,{'CPM','RAD'})} = 'CPM'
    option (1,:) char {mustBeMember(option,{'TRIM','ALL','SameArea'})} = 'ALL'
    pl (1,1) logical = true
end

switch(source)
    
    case 'CPM'
        ENSEMBLENO = getEnsNos();
        rmInd = CELLALL.Area<=1;
        CountAll = [];
        for ensNo = 1:length(ENSEMBLENO)
            C0 = getCounts(CELLALL,rmInd,ensNo,E);
            % C0 = conv2(C0,ones(3,3)/9,'same');
            CountAll(:,:,ensNo) = reshape(C0,[size(C0),1]);
        end
        Count = squeeze(nanmean(CountAll,3));
        
    case 'RAD'
        % XEDGES,YEDGES
        EEdges = 0.5:1:size(E,1)+0.5;
        NEdges = 0.5:1:size(N,2)+0.5;
        rmInd = CELLALL.Area<=1;
        Count = histcounts2(CELLALL.Centroid(~rmInd,2),CELLALL.Centroid(~rmInd,1),EEdges,NEdges);
        Count = Count/12;
        % Count = conv2(Count,ones(3,3)/9,'same');
        CountAll = Count;
end

UKMap = getUKMap();

if strcmpi(option,'TRIM')
    
    switch(source)
        case 'CPM'
            in = getTrimTag('unit','km','product','cpm2.2');
        case 'RAD'
            in = getTrimTag('unit','km','product','radar2.2');
    end
    Count(~in) = NaN;
    
elseif strcmpi(option,'ALL')
    
elseif strcmpi(option,'SameArea')
    
    switch(source)
        case 'CPM'
            in = getTrimTag('unit','km','product','cpm2.2');
        case 'RAD'
            in = getTrimTag('unit','km','product','radar2.2');
    end
    [E0,N0,Count] = remapCount(E,N,Count,in);
    
    if numel(size(CountAll))==3
        CountAll2 = [];
        for ensNo = 1:size(CountAll,3)
            [~,~,D] = remapCount(E,N,squeeze(CountAll(:,:,ensNo)),in);
            CountAll2(:,:,ensNo) = reshape(D,[size(D),1]);
        end
        CountAll = CountAll2;
    else
        CountAll = Count;
    end
    E = E0;
    N = N0;
else
    
end


f = @(mat)mat;%conv2(mat,ones(3,3)/9,'same');% imresize(mat,0.5,'box');

if pl
    showVal = f(Count);
    showVal(showVal == 0) = NaN;
    
    pcolor(f(E),f(N),showVal);
    hold on;
    shading flat
    
    upperl = 2;%prctile(Count(Count>0),90);
    caxis(gca,[0,upperl])
    cptcmap('OrRd_09','mapping', 'scaled','ncol',9,'flip',false);
    alpha(0.8);
    xlabel('Easting/Km');
    plot(UKMap.borderE/1000,UKMap.borderN/1000,'k-');
    ylabel('Northing/Km')
    XYWH = [150,0,420,490];
    set(gcf,'units','points','position',XYWH);
    
    if numel(seasons) == 1
        text(-50,1300,sprintf('Cell Number-%s(%s)',getSeasonName(seasons),source),...
            'fontweight','bold','fontsize',14,...
            'verticalalignment','top',...
            'horizontalalignment','left');
    else
        text(-50,1300,sprintf('Cell Number  (%s)',source),...
            'fontweight','bold','fontsize',14,...
            'verticalalignment','top',...
            'horizontalalignment','left');
    end
    % colorbar;
    xlim([-50,750]);
    ylim([0,1300])
    axis off
end

    function C0 = getCounts(CELLALL,rmInd,ensNo,E)
        C0 = histcounts2(CELLALL.Centroid(~rmInd & CELLALL.ensNo==ensNo,2),...
            CELLALL.Centroid(~rmInd & CELLALL.ensNo==ensNo,1),size(E));
        C0 = C0/20;
    end
    function [E,N,Count] = remapCount(E,N,Count,in)
        
        Count(~in) = NaN;
        
        dx = 2.2;
        [E_same,N_same] = meshgrid(-280:dx:817.8, -226:dx:1179.8);
        
        Count0 = griddata(E(:),N(:),Count(:),E_same(:),N_same(:),'nearest');
        Count0 = reshape(Count0,size(E_same));
        Count = Count0;
        E = E_same;
        N = N_same;
    end

end









