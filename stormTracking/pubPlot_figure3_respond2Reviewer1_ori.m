% ----------------------------------------------------------------------- %
% THIS SCRIPT IS TO VISUALIZE SEVERAL STATISTICS FOR CONVECTIVE STROM DURING
% SUMMERTIME
%
% ... # some necessary description # ...
%
% # Spatial Characteristics
% Here storms with PMax were identified, and then density center were
% identified for plotting
%
% # Smoothing
% Convolution was used in several steps to locate the Pmax. But of course
% on the final images showing the averaged pattern, we are using the raw
% simulation output from CPM.
%
%
%
% @ Yuting Chen
% Imperial College London
% yuting.chen17@imperial.ac.uk
% ----------------------------------------------------------------------- %
clear;clc;

% Several Config
close all
setFigureProperty('Paper');
global regionName savePath region ensNo
ENSEMBLENO = getEnsNos();
durThre = [];% 6;

Num = 40;%
for regionName = {'CPM_NW','CPM_NE','CPM_S'}% {'SCO','WAL','EUK'}% 

    regionName = regionName{1};
    region = getfield(REGIONS_info(),regionName);
    
    x_name = 'rpmax';%'rspeed';%'rsize';%
    savePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\StormTracking';
    %
    PERIODS = {'1980-2000','2060-2080'}; % {'2007-2018'};% ;'2020-2040',
    RAs = [];
    ha_PaFu = tight_subplot(1,3,[0.01 0.01],[.19 0.01],[.10 .15]);
    % Get data (opt: plot)
    pi = 1;
    summaryTable = [];
    figure(1)
    for Period = PERIODS
        
        Period = Period{1};
%         %%
%         Ra = [];
%         for ensNo = ENSEMBLENO([1:12])
%             ensNo = ensNo{1};
%             [STATS0,T0] = get4Plot(Period,ensNo,Num,durThre);
%             %%
%             Ra0 = plot_40CS(STATS0,T0,Period,durThre,ensNo);
%             Ra = cat(3,Ra,Ra0);
%             close all
%             
%             % plot_averageRa(Ra0,Period,{ensNo},durThre);
%             close all
%         end
%         Config = getConfig(upper(regionName),6,Period,ensNo);
%         save(sprintf('%s%sCS_%s_%s_4PubPlot_%03dStorms.mat',Config.saveIt.path,...
%             filesep,regionName,Period,Num*12),'Ra');
        
        
        %%
        axes(ha_PaFu(pi))
        ensNo = '15';
        Config = getConfig(upper(regionName),6,Period,ensNo);
        load(sprintf('%s%sCS_%s_%s_4PubPlot_%03dStorms.mat',Config.saveIt.path,...
            filesep,regionName,Period,Num*12),'Ra');
        
        [summaryTable{pi}] = plot_averageRa(Ra,Period,ENSEMBLENO,durThre);
        
        RAs{pi} = Ra;
        pi = pi+1;
    end
    
    c = colorbar('Location','SouthOutside','Position',[0.1 0.055 0.5 0.05]);
    c.Ticks = [1,20,40,60,80];c.TickLength = 0.02;
    c.TickLabels = strcat(c.TickLabels,'mm/h');
    
    ax = gca;
    ylabel('');
    ax.YTickLabel = [];
  
    figure(3)
    ratio_cov = summaryTable{2}.spatialCoverage./summaryTable{1}.spatialCoverage;
    ratio_vol = summaryTable{2}.totalVolume./summaryTable{1}.totalVolume;
    ratio_pi = summaryTable{2}.peakIntensity./summaryTable{1}.peakIntensity;
    % plot(ratio_pi,ratio_vol,'o');hold on;
    
    %     y = ratio_cov;
    xrange = summaryTable{2}.peakIntensity_range./summaryTable{1}.peakIntensity_range;
    yrange = summaryTable{2}.totalVolume_range./summaryTable{1}.totalVolume_range;
    x = nanmedian(xrange);
    y = nanmedian(yrange);
    %     yrange = summaryTable{2}.spatialCoverage_range./summaryTable{1}.spatialCoverage_range;
    yneg = std(yrange);
    ypos = std(yrange);
    xneg = std(xrange);
    xpos = std(xrange);
    errorbar(x,y,yneg,ypos,xneg,xpos,'o');hold on
    

    xlabel('P_{max}^{fut}/P_{max}^{pas}-1')
    ylabel('P_{mean}^{fut}/P_{mean}^{pas}-1')
    legend boxoff
    axis([0.5,2,0.5,2]);
    
    axis('square')
    % close all
    figure(2)
    
end
figure(3);
hold on;
% plot([0.5,2],[0.5,2],'k--','linewidth',0.5);% 
plot([1,1],[0.5,2],'k--','linewidth',0.5);
plot([0.5,2],[1,1],'k--','linewidth',0.5);
legend('UK-NW','UK-NE','UK-S','Location','Southeast');

% AUXILLARY FUNCTION
function Col = getColor(period)
switch(period)
    case '1980-2000'
        Col = [66 146 199]/255;
    case '2060-2080'
        Col = [240  60  43]/255;
end
end
function hand=plotOneIntArea(DAT,RD,RA,colo);
% DAT0 = DAT(:);
% RD0 = floor(RD(:));
% T = table(DAT0,RD0);

RD0 = repmat(reshape(RD,[size(RD),1]),[1,1,size(RA,3)]);
RD0 = RD0(:);
DAT0 = RA(:);
T = table(DAT0,RD0);
A = grpstats(T,'RD0',{'median',@(x)quantile(x,.75),@(x)quantile(x,.25)});
xi = A.RD0;

try
% % #fit version#
% xx = xi;
% yy = A.median_DAT0;
% ft = fittype('pm*exp(-x/k)+d');
% Mfit = fit(xx,yy,ft,'StartPoint',[nanmax(yy)+0.1,9,0.06]);
% k = Mfit.k;
% d = Mfit.d;
% pmax = Mfit.pm;
% ci = confint(Mfit,0.95);
% rr = linspace(0.1,nanmax(xi),200);
% xx = pmax.*exp(-rr./k)+d;
% yy = 2*pi*(rr).^2;
% hand = plot(xx,yy,'-','color',colo);

% #no-fit version#
% xx = A.mean_DAT0;
% yy = 2*pi*(xi+0.5).^2;
hand = plot(A.median_DAT0,0.1+2*pi*(xi).^2,'--','color',colo);
catch me
    % #no-fit version#
    % xx = A.mean_DAT0;
    % yy = 2*pi*(xi+0.5).^2;
    hand = plot(A.median_DAT0,0.1+2*pi*(xi).^2,'--','color',colo);
end
hsimran = fill([A.Fun3_DAT0;flip(A.Fun2_DAT0)],0.1+2*pi.*([xi;flip(xi)]).^2,colo,...
    'LineStyle','none');
alpha(0.3)
end

function Ra = plot_40CS(STATS0,T0,Period,durThre,ensNo)

global regionName savePath region
figure;
Ra = zeros(51,51,0);
ha = tight_subplot(5,ceil(length(STATS0)/5),[0.0 0.0],[.05 0.05],[.05 .15]);
for i = 1:length(STATS0)
    R = STATS0{i};
    thisTime = T0{i};
    
    % # system centre #
    % # opt 1 #
    % R_temp = conv2(R,ones(3,3)/9,'same');
    % [maxi,maxj] = find(R_temp==nanmax(R_temp(:)),1);
    % fprintf('%3dmm/h\n',nanmax(R_temp(:)));
    
    % # opt 2 # #seems better#
    % This one gives us the weighted centroid for area having highest Pmax.

    R_temp = R*0; 
    R_temp([1+25:size(R,1)-25],[1+25:size(R,2)-25]) = R([1+25:size(R,1)-25],[1+25:size(R,2)-25]);

    % pick up heaviest point.
    surVec = [1,1,1;1,0,1;1,1,1];
    R_temp(R_temp<5) = 0;
    surB5 = conv2(R_temp,surVec,'same');% have at least one surronding bigger than 5mm/h
    R_temp(surB5<1) = 0; % ensure at least 2 grids > 5mm/h
    allVec = [1,1,1;1,2,1;1,1,1];%%%%%%%%%%%%%%%0515Final
    R_temp = conv2(R_temp,allVec/10,'same');%%%%%%%%%%%%%%%%%%%%%%0515Final
    [maxi,maxj] = find(R_temp==nanmax(R_temp(:)),1);
    % R_temp(R_temp<floor(nanmax(R_temp(:)))) = 0;
    % stats = regionprops('table',logical(R_temp),(R_temp),'WeightedCentroid','Area','MaxIntensity');
    
    
    % system centroid
    % R_temp(R_temp<5) = 0;
    % stats = regionprops('table',logical(R_temp),(R_temp),'WeightedCentroid','Area','MaxIntensity');
    % stats(stats.Area==1,:) = [];
    % maxi = round(stats.WeightedCentroid(stats.MaxIntensity == max(stats.MaxIntensity),2));
    % maxj = round(stats.WeightedCentroid(stats.MaxIntensity == max(stats.MaxIntensity),1));
    R(R==0) = NaN;

    R = squeeze(R([max(maxi-25,1):min(size(R,1),maxi+25)],...
        [max(1,maxj-25):min(size(R,2),maxj+25)]));
    
    if ~any(~(size(R)==[51,51]))
        R(isnan(R))=0;[maxi,maxj] = deal(26);
        Ra = cat(3,Ra,R);
    else
        R = NaN(51,51);
        Ra = cat(3,Ra,R);
%         % # In some rare case #.
%         % centroid point might be very close to boundary of area.
%         % In this case, current setting is to enlarge area and then search
%         % for the complete image pattern again.
%         elNo = 20;
%         if strcmp(region.Name,'bigEUK')
%             elNo = 5;
%         end
%         imageNo = (thisTime.Day-1)*24+thisTime.Hour+1;
%         [A,~,~] = readCPM_nc_1_wider(region,thisTime.Year,thisTime.Month,ensNo,imageNo,elNo);
%         A(A<1/32)=0;
%         R = A;
%         maxi = maxi+elNo;
%         maxj = maxj+elNo;
%         R(R==0) = NaN;
%         R = squeeze(R([max(maxi-25,1):min(size(R,1),maxi+25)],...
%             [max(1,maxj-25):min(size(R,2),maxj+25)]));
%         % R = R*100;
%         if any(~(size(R)==[51,51]))
%             1;
%             continue;
%         else
%             R(isnan(R))=0;
%             Ra = cat(3,Ra,R);
%         end
    end
    if ~any(~(size(R)==[51,51]))
        axes(ha(i));hold on;
        R(R==0)=NaN; pcolor(1:51,1:51,R);
        % [maxi,maxj] = find(R == nanmax(R(:)));
        maxj = 26;maxi = 26;
        % plot(maxj,maxi,'r+','linewidth',2);
        % xlim(maxj+[-50,50]);
        % ylim(maxi+[-50,50]);
        shading flat;box on;
        cptcmap('cw1-013','mapping','scaled','ncol',100,'flip',true);
        caxis([0,50]);
        format short g
        thisTime.Format='dd-MMM-uuuu HH:mm';
        text(25,7,[string(thisTime),sprintf('%smm/h',num2str(round(R(26,26))))],...
            'horizontalAlignment','center','fontsize',8,'fontweight','bold');
        text(1,50,sprintf('(%d)',i),...
            'horizontalAlignment','left',...
            'verticalAlignment','top',...
            'fontsize',8,'fontweight','bold');
        axis([1,51,1,51])
        ax = gca;ax.XTick = [];ax.YTick = [];
        axis('square')
    end
end

% c = colorbar(gca);
c = colorbar('location','Manual', 'position', [0.92 0.1 0.015 0.8]);
c.Ticks = 0:10:50;
c.TickLabels = strcat(c.TickLabels,'mm/h');
saveName = sprintf('%s-40CSs-highPMax_%s_%s-Thre%02d',regionName,Period,ensNo,durThre);
savePlot([savePath,filesep,saveName],'XYWH',[150,0,909,505],'needreply','N','onlyPng',true);
end

function [hand,xi,yi,Mfit] = plotReduction(DAT,RD,RA,colo);

[Mfit] = deal(NaN);
imN = 5;
RA = imresize3(RA,'scale',[imN,imN,1],'method','box');
x0 = 2.2*linspace(-25,25,51*imN);
y0 = 2.2*linspace(-25,25,51*imN);
[XX,YY] = meshgrid(x0,y0);
RD = sqrt(XX.^2+YY.^2);

RA0 = reshape(RA,[size(RA,[1,2]),size(RA,3)/12,12]);
RD0 = repmat(reshape(RD,[size(RD),1]),[1,1,size(RA0,[3,4])]);
ensNo = repmat(reshape(1:12,[1,1,1,12]),[size(RA0,[1:3]),1]);
RD0 = round(RD0(:));
DAT0 = RA0(:);
ensNo = ensNo(:);
T = table(DAT0,RD0,ensNo);

A = grpstats(T,{'RD0','ensNo'},{'mean',@(x)quantile(x,.75),@(x)quantile(x,.25)});
RD0 = A.RD0;
DAT0 = A.mean_DAT0;
ensNo = A.ensNo;

T = table(DAT0,RD0,ensNo);
A = grpstats(T,{'RD0'},{'median',@(x)quantile(x,1),@(x)quantile(x,0)});
xi = A.RD0;
yi = A.median_DAT0./nanmax(A.median_DAT0);
hand = plot(pi*xi.^2,yi,'-','color',colo);
hsimran = fill([pi*xi.^2;flip(pi*xi.^2)],[A.Fun3_DAT0;flip(A.Fun2_DAT0)],colo,...
    'LineStyle','none');

xx = xi;
yy = A.median_DAT0;
ft = fittype('pm*exp(-x/k)+d');
Mfit = fit(xx,yy,ft,'StartPoint',[nanmax(yy)+0.1,9,0.06],'Exclude', xx > 100);
k = Mfit.k;
d = Mfit.d;
pmax = Mfit.pm;
ci = confint(Mfit,0.95);
xx = linspace(nanmin(xx),nanmax(xx),100);
yy = pmax.*exp(-xx./k)+d;
hold on
hand = plot(pi*xx.^2,yy,'--','color',colo);

alpha(0.3)
end


function [hand,handran,xi,yi,Mfit,T] = plotOneRadialP(DAT,RD,RA,colo);

% DAT0 = DAT(:);
% RD0 = floor(RD(:));
% T = table(DAT0,RD0);
[Mfit] = deal(NaN);
imN = 5;
RA = imresize3(RA,'scale',[imN,imN,1],'method','box');
x0 = 2.2*linspace(-25.5,25.5,size(DAT,1)*imN);
y0 = 2.2*linspace(-25.5,25.5,size(DAT,1)*imN);
[XX,YY] = meshgrid(x0,y0);
RD = sqrt(XX.^2+YY.^2);

RA0 = reshape(RA,[size(RA,[1,2]),size(RA,3)/12,12]);
RD0 = repmat(reshape(RD,[size(RD),1]),[1,1,size(RA0,[3,4])]);
ensNo = repmat(reshape(1:12,[1,1,1,12]),[size(RA0,[1:3]),1]);
RD0 = round(RD0(:)/2.2)*2.2;
RD0 = (RD0(:));
DAT0 = RA0(:);
ensNo = ensNo(:);
T = table(DAT0,RD0,ensNo);

A = grpstats(T,{'RD0','ensNo'},{'mean',@(x)quantile(x,.75),@(x)quantile(x,.25)});
RD0 = A.RD0;
DAT0 = A.mean_DAT0;
ensNo = A.ensNo;

T = table(DAT0,RD0,ensNo);
A = grpstats(T,{'RD0'},{'median',@(x)quantile(x,1),@(x)quantile(x,0)});
xi = A.RD0;
yi = A.median_DAT0;
hand = plot(xi,yi,'-','color',colo);hold on;
handran = fill([xi;flip(xi)],[A.Fun3_DAT0;flip(A.Fun2_DAT0)],colo,...
    'LineStyle','none');
try
xx = xi;
yy = A.median_DAT0;
ft = fittype('pm*exp(-x/k)+d');
Mfit = fit(xx,yy,ft,'StartPoint',[nanmax(yy)+0.1,9,0.06],'Exclude', xx > 100);
k = Mfit.k;
d = Mfit.d;
pmax = Mfit.pm;
ci = confint(Mfit,0.95);
xx = linspace(nanmin(xx),nanmax(xx),100);
yy = pmax.*exp(-xx./k)+d;
hold on
% hand = plot(xx,yy,'--','color',colo);
%     hand = plot(xi,yi,'-','color',colo);
catch me
%     hand = plot(xi,yi,'-','color',colo);
end
% hsimran = fill([xi;flip(xi)],[A.Fun3_DAT0;flip(A.Fun2_DAT0)],colo,...
%     'LineStyle','none');
alpha(0.3)
end

function [summary] = plot_averageRa(Ra,Period,ENSEMBLENO,durThre)
global regionName savePath ensNo region
fprintf('%02d snapshots were used\n',size(Ra,3));
% figure;
numRa = size(Ra,3);
summary = struct;

Ra = reshape(Ra,[size(Ra,[1,2]),size(Ra,3)/length(ENSEMBLENO),length(ENSEMBLENO)]);
summary.spatialCoverage_range = squeeze(nansum(squeeze(nanmean(Ra,3))>5,[1,2])*2.2*2.2);
summary.totalVolume_range = squeeze(nansum(squeeze(nanmean(Ra,3)),[1,2])*10^-3/3600*2200*2200);
summary.peakIntensity_range = squeeze(nanmax(squeeze(nanmean(Ra,3)),[],[1,2]));

Ra = squeeze(nanmedian(nanmean(Ra,3),4));
summary.totalVolume = nansum(Ra(:))*10^-3/3600*2200*2200;
summary.spatialCoverage = nansum(Ra(:)>5)*2.2*2.2;
summary.peakIntensity = nanmax(Ra(:));

Ra = Ra(11:41,11:41);
scaleN = 1;%22;
Ra = -0.1+exp(imresize(log(Ra+0.1),scaleN,'lanczos2'));
% Ra = -0.1+exp(imresize(log(Ra+0.1),scaleN,'bilinear'));
x = ((1:size(Ra,1))-(size(Ra,1)+1)/2)*2.2/scaleN;
y = ((1:size(Ra,2))-(size(Ra,1)+1)/2)*2.2/scaleN;
Ra_temp = Ra;

Ra_temp(Ra_temp<0.1) = NaN;%%%%%%%% 
pcolor(x,y,Ra_temp);hold on;
grid off
[C,h] = contour(x,y,Ra,[5,10,20],'Color',[0.1 0.1 0.1],'Linewidth',0.5,'Linestyle','-');
clabel(C,h,'FontSize',10)%,'Color','red')
shading flat
cptcmap('cw1-013','mapping','scaled','ncol',32,'flip',true);
caxis([0,80])
grid on;
%
% c = colorbar;
% c.Ticks = [1,20,40,60,80];c.TickLength = 0.02;
% c.TickLabels = strcat(c.TickLabels,'mm/h');
%
% plot(size(Ra,1)/2+1,size(Ra,2)/2+1,'r+')
% set(gca,'TickDir','out');

xlabel('Distance from Centre /km');
ylabel('Distance from Centre /km');
ax = gca;hold on;
plot(mean(ax.XLim)*ones(2,1)+1.1,ax.YLim+1.1,'k:','linewidth',0.5);
plot(ax.XLim+1.1,mean(ax.YLim)*ones(2,1)+1.1,'k:','linewidth',0.5);
text(ax.XLim(1),ax.YLim(2),sprintf('%s\n%s',...
    getLabel(regionName),Period),...
    'horizontalalignment','left',...
    'verticalalignment','top','fontsize',10,'Interpreter','latex');
text(ax.XLim(2),ax.YLim(2),sprintf('%4.1f$m^3 s^{-1}$',...
    summary.totalVolume),...
    'horizontalalignment','right',...
    'verticalalignment','top','fontsize',10,'Interpreter','latex');
% plot(20*ones(2,1),ax.YLim,'k:','Color',[0.1 0.1 0.1],'linewidth',0.5);
% plot(-20*ones(2,1),ax.YLim,'k:','Color',[0.1 0.1 0.1],'linewidth',0.5);
% plot(ax.XLim,-20*ones(2,1),'k:','Color',[0.1 0.1 0.1],'linewidth',0.5);
% plot(ax.XLim,20*ones(2,1),'k:','Color',[0.1 0.1 0.1],'linewidth',0.5);
ax.YTick = [-50,-20,0,20,50];
ax.XTick = [-50,-20,0,20,50];
hold off;
box on
axis('square')

% if iscell(ENSEMBLENO) && length(ENSEMBLENO)>1
%     saveName = sprintf('%s-%03dCSs-averaged_%s_%s-%s-Thre%03d',regionName,numRa,Period,...
%         ENSEMBLENO{1},ENSEMBLENO{end},durThre);
% else
%     saveName = sprintf('%s-%03dCSs-averaged_%s_%s-Thre%03d',regionName,numRa,Period,...
%         ensNo,durThre);
% end
% savePlot([savePath,filesep,saveName],'XYWH',[150,0,350,250],'needreply','N',...
%     'onlyPng',true);
end


function [Rain,scaleF,region] = readCPM_nc_1_wider(region,year,mon,ensNo,imageNo,elNo);

fileName = sprintf('pr_rcp85_land-cpm_uk_2.2km_%s_1hr_%04d%02d01-%04d%02d30.nc',...
    ensNo,year,mon,year,mon);
try
    filePath = ['K:/UkCp18/',ensNo,'/'];
    listaRain = [filePath,fileName];
    A = ncinfo(listaRain);
catch
    filePath = ['K:/UkCp18_FutureTemp/',ensNo,'/'];
    listaRain = [filePath,fileName];
    A = ncinfo(listaRain);
end

LAT=ncread(listaRain,'latitude');
LON=ncread(listaRain,'longitude');
[E, N] = ll2os(LAT, LON);
E = E/1000;
N = N/1000;
[region.i,region.j] = getRegionIJ(E,N,region.minE,region.minN);

Rain = squeeze(ncread(listaRain,'pr',...
    [region.i-elNo,region.j-elNo  imageNo,  1],...
    [region.dimE+elNo*2,region.dimN+elNo*2    1,  1]));

scaleF = 1;
end

function str = getLabel(regionName)
switch(regionName)
    case 'SCO'
        str = 'NUK';
    case 'EUK'
        str = 'SEUK';
    case 'WAL'
        str = 'SWUK';
    case 'CPM_NW'
        str = 'NW-UK';
    case 'CPM_NE'
        str = 'NE-UK';
    case 'CPM_S'
        str = 'S-UK';
end
end
function [A,T] = get4Plot(Period,ensNo,Num,durThre)

global regionName region
[RE,TE] = deal([]);
% for MON = 6:8
%     Config = getConfig(upper(regionName),MON,Period,ensNo);
%     Dat = load([Config.saveIt.path,filesep,sprintf('CS_%s_%s_FIELD_%02d_%s.mat',regionName,Period,Config.Month,ensNo)],...
%         'RE','TE','Config');
%     RE = cat(1,RE,Dat.RE);
%     TE = cat(1,TE,Dat.TE);
%     % STATS = [STATS;Dat.RE];
% end

Config = getConfig(upper(regionName),6,Period,ensNo);
Dat = load([Config.saveIt.path,filesep,sprintf('CS_%s_%s_FIELD_%02d-%02d_%s.mat',regionName,Period,6,8,ensNo)],...
    'RE','TE','Config');
RE = Dat.RE;
TE = Dat.TE;
len = cell2mat(cellfun(@(a)size(a,3),RE,'UniformOutput', false));
if ~isempty(durThre)
RE = RE(len>durThre);
TE = TE(len>durThre);
end
% STATS = [STATS;Dat.RE];

% trim out non-land area for each event
[E,N] = getEN(region);
UKMap = load('H:\CODE_MATLAB\SpatialTemporalDATA\shapeFileFolder\UKBorderGrid.mat');
[in] = inpolygon(E,N,UKMap.borderE/1000,UKMap.borderN/1000);
RE = cellfun(@(a)trimIt(a,in),RE,'UniformOutput', false);

% only search for the max value for central part
STATS_rpmax = cell2mat(cellfun(@(a)nanmax(reshape(a(1+25:size(RE{1},1)-25,1+25:size(RE{1},2)-25,:),...
    1,[])),...
    RE,'UniformOutput', false));
[~,I] = sort(STATS_rpmax,'descend');
A = RE(I(1:Num));
STATS_rpmax(I(1:Num));
T = TE(I(1:Num));

[A,T] = cellfun(@(a,t)findRep(a,t,[1+25:size(RE{1},1)-25],[1+25:size(RE{1},2)-25]),A,T,'UniformOutput', false);
    function a = trimIt(a,in)
        a(repmat(~in,[1,1,size(a,3)])) = 0;%NaN;
    end
    function [A,T] = findRep(As,Ts,Rii,Rjj)
        % # Opt 1 #
        % We will take the one with biggest coverage of wet as the
        % representaed images of this events.
        % war = nanmean(reshape(As,[prod(size(As,[1,2])),size(As,3)])>1,1);
        % maxi = find(war == nanmax(war));
        % A = squeeze(As(:,:,maxi));
        % T = Ts(maxi);
        %
        % # Opt 2 #
        % Here a 3*3 convolution was used to smooth spatial
        % field, and then one image having highest Pmax was picked up.
        
        % # Opt 3 #
        As0 = As(Rii,Rjj,:);
        
        surVec = [1,1,1;1,0,1;1,1,1];
        allVec = [1,1,1;1,2,1;1,1,1];
        for i = 1:size(As0,3)
            atemp = squeeze(As0(:,:,i));% 
            % atemp = conv2(squeeze(As0(:,:,i)),ones(3,3)/9,'same');
            % As0(:,:,i) = atemp;
            
            % exclude the case: only one single point as pmax but all
            % surrondings are smaller than 5 mm/h.
            atemp(atemp<5) = 0;
            surB5 = conv2(atemp,surVec,'same');% have at least one surronding bigger than 5mm/h
            atemp(surB5<1) = 0;% ensure at least 2 connected is higher than 5mm/h
            
            atemp = conv2(atemp,allVec/10,'same');%%%%%%%%%%0515Final
            As0(:,:,i) = atemp;
            
        end
        % As0 = As;
        pm = squeeze(nanmax(nanmax(As0,[],1),[],2));
        indmax = find(pm==nanmax(pm(:)),1);
        A = squeeze(As(:,:,indmax));
        T = Ts(indmax);
        if ~isRegularRegion(region)
            angle = atan((region.E(2)-region.E(1))/(region.N(2)-region.N(1)))*180/pi;
            A = imrotate(A,-angle,'nearest','loose');%rotate it back
        end
    end
end
