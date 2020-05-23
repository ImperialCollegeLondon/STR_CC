% ----------------------------------------------------------------------- %
% THIS FILE IS TO INVESTIGATE IF CLIMATE MODELS GAVE A GOOD CORRELATION 
% BETWEEN TOPOGRAPHY AND RAINFALL GIVEN A DEFINED CROSS-SECTION
% 
% DATA USED:
%           Climate Data: UKCP18 CPM 2.2
%           RainGauge Data: Monthly GEAR from CEH
%           Radar Data: MetOffice NIMROD Radar Composite 5min
%           Topo Data: Ordance Survey TERRAIN50 DTM
% 
% @ Yuting CHEN
% yuting.chen17@imperial.ac.uk
% Update: 2020.01.27
% ----------------------------------------------------------------------- %



clear;clc

REGIONS = REGIONS_info();
region = REGIONS.Scotland;% Westuk; % wales



%% Result from GEAR-Month

[GEAR,X_coor,Y_coor,DIST] = getRegionalGEAR(region,'month');

%%

UKMap = load('H:\CODE_MATLAB\SpatialTemporalDATA\shapeFileFolder\UKBorderGrid.mat');
in = inpolygon(X_coor,Y_coor,UKMap.borderE/1000,UKMap.borderN/1000);
GEAR(:,~in) = NaN; 

%%
[gearIJ] = find(DIST>1000);
GEAR(:,gearIJ) = NaN;

%% Result from RADAR

[RAD,E_rad,N_rad] = getMonthRadar(region);


%% Result from CPM

[CPM,E,N] = getMonthCPM(region);% {mon,ensNo}[E,N,yr]


%% Result for TOPO

[Topo,E,N] = getTopo(region);


%% DEFINE CROSS-SECTION
p0 = [206700,735400];
p1 = [314500,627600];
dd = linspace(0,1,50*3/2.2);
pCross = transpose(p0' + dd.*(p1'-p0'));


%% plot TOPO and Cross-section

pcolor(E,N,Topo);shading flat
cptcmap('DEM_print', 'mapping','direct');
hold on;
plot(pCross(:,1),pCross(:,2),'r','linewidth',10)
colorbar

filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\Scotland_OneCrossSection';
filename = [filePath,filesep,'Topo_Location'];
savePlot(filename,'XYWH',[150,0,800,700],'NeedReply','Y');



%% get CROSS_TOPO

F = scatteredInterpolant(E(:),N(:),Topo(:),'nearest');
vq = F(pCross(:,1),pCross(:,2));


%% filter out topo outside uk boundary

UKMap = load('H:\CODE_MATLAB\SpatialTemporalDATA\shapeFileFolder\UKBorderGrid.mat');
in = inpolygon(E_rad*1000,N_rad*1000,UKMap.borderE,UKMap.borderN);
Topoco = imresize(Topo,size(E_rad),'nearest');
Topoco(~in) = NaN;

%% PLOT TOPO-RAIN


UKMap = load('H:\CODE_MATLAB\SpatialTemporalDATA\shapeFileFolder\UKBorderGrid.mat');
in_rad = inpolygon(E_rad*1000,N_rad*1000,UKMap.borderE,UKMap.borderN);

mon = 2;
ha = tight_subplot(2,3,[.1 .1],[.15 .05],[.1 .1]);

for enNo = 1:12
    
    CPM{mon,enNo}(~in_rad) = NaN;
    RAD(mon,~in_rad) = NaN;
    
    plotTOPO_RAIN(ha,...
        Topoco,CPM{mon,enNo},...
        Topoco,squeeze(RAD(mon,:,:)),...
        imresize(Topoco,size(squeeze(GEAR(mon,:,:))),'nearest'),...
        transpose(squeeze(GEAR(mon,:,:)/eomday(1990,mon)/24)),enNo == 12);
end

mon = 8;
for enNo = 1:12
    
    CPM{mon,enNo}(~in_rad) = NaN;
    RAD(mon,~in_rad) = NaN;
    
    plotTOPO_RAIN(ha(4:6),...
        Topoco,CPM{mon,enNo},...
        Topoco,squeeze(RAD(mon,:,:)),...
        imresize(Topoco,size(squeeze(GEAR(mon,:,:))),'nearest'),...
        transpose(squeeze(GEAR(mon,:,:)/eomday(1990,mon)/24)),enNo == 12);
end


arrayfun(@(han)xlabel(han,''),ha(1:3));
arrayfun(@(han)grid(han,'minor'),ha);

filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP';
filename = [filePath,filesep,'Topo_Scotland_wholeBox'];
savePlot(filename,'XYWH',[150,0,700,400]);




%% AUXILLARY FUNCTION


function plotTOPO_RAIN(ha,TopoM,CPMM,...
    TopoM2,RADM,...
    TopoM3,GEARM,tag)

arguments
    ha (:,1)
    TopoM (:,:) double = NaN(50,50);
    CPMM (:,:) double = NaN(50,50);
    TopoM2 (:,:) double = NaN(50,50);
    RADM (:,:) double = NaN(50,50);
    TopoM3 (:,:) double = NaN(50,50);
    GEARM (:,:) double = NaN(50,50);
    tag (1,1) logical = 0;
end

cmap = gray(6);
cmap(1,:) = [];
setFigureProperty('Subplot3')

axes(ha(1))
plotG(TopoM,CPMM,1,'CPM2.2');
alpha(0.1)
hold on;
grid minor

if tag
axes(ha(2))
plotG(TopoM2,RADM,2,'RAD');
hold on;
grid minor

axes(ha(3))
plotG(TopoM3,GEARM,3,'GEAR');
hold on;
grid minor
end
    function plotG(xMat,yMat,cmapNo,product)
        x = xMat(:);
        y = reshape(yMat*24,[],1);
        ind = (~isnan(x)) & (~isnan(y));
        x = x(ind);
        y = y(ind);
        
        fitresult = fit(x,y,'poly1');
        p22 = predint(fitresult,unique(sort(x)),0.95,'functional','on');
        
        scatter(x,y,10,cmap(cmapNo,:),'filled');
        hold on;
        
        plot(fitresult,'k-')
        hold on;
        plot(unique(sort(x)),p22,'k--','Linewidth',1);
        hold on;
        
        box on;
        xlabel('DTM /m')
        ylabel('Rainfall mm/h')
        ylim([0.1,0.5]*24)
        
        XLIM = xlim;
        YLIM = ylim;
        if tag
        text(XLIM(2)*0.95,YLIM(2)*0.95,product,'HorizontalAlignment','right',...
            'VerticalAlignment','top','Backgroundcolor','w');
        end
        alpha(0.3)
        legend off
    end

end


function [Topo,E,N] = getTopo(region)

TERRAIN = load(['K:\UK_shape\DTM50.mat'],'DTM50','Eno','Nno');

Nind_max = findClose(region.minN*1000,TERRAIN.Nno);
Nind_min = findClose(1000*(region.minN+region.dx*(region.dimN-1)),TERRAIN.Nno);

Eind_min = findClose(region.minE*1000,TERRAIN.Eno);
Eind_max = findClose(1000*(region.minE+region.dx*(region.dimE-1)),TERRAIN.Eno);

Topo = TERRAIN.DTM50(Nind_max:-1:Nind_min,...
    Eind_min:Eind_max);

[NNno,EEno] = meshgrid(TERRAIN.Nno(Nind_max:-1:Nind_min),TERRAIN.Eno(Eind_min:Eind_max));
E = EEno;
N = NNno;
Topo = transpose(Topo);

% % regrid to 2.2km reso.
% E = region.minE:2.2:region.minE+(region.dimE-1)*2.2;
% N = region.minN:2.2:region.minN+(region.dimN-1)*2.2;
% [N,E] = meshgrid(N*1000,E*1000);% not same as X-coor,Y-coor
% tic
% Topo = griddata(EEno(:),NNno(:),Topo(:),E,N,'natural');
% toc



    function indX = findClose(x,Xvec)
        dis = abs(Xvec - x);
        indX = find(dis == min(dis(:)),1);
    end

end




