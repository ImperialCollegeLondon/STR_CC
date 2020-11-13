% This script is to plot IDAF based on annual maxima ananlysis
%
% Data: 1980-2000 UKCP CPM rainfall or 2060-2080 UKCP CPM rainfall
%
% ref: ...2002
% based on C. De Michele 2011 func.6 in Journal of Hydrology
%
% @ Yuting CHEN
% yuting.chen17@imperial.ac.uk
% update: 2020.08.06


clear;clc
close all
warning on

futPeriod = '2060-2080';
pasPeriod = '1980-2000';
REGIONS = {'CPM_NW','CPM_NE','CPM_S'}; % {'SCO','WAL','EUK'};
dataInfo = getDataInfo(REGIONS{1},pasPeriod);
rt = exp(linspace(log(2),log(40),40));
area = 2.^(unique(round(5*log2(1:8000)))/5);
figForm = struct; % format of output figure.
[figForm.rt,figForm.area] = meshgrid(rt, area);
itag = 1; Zs = [];


% ha = tight_subplot(1,3,[.05 .05],[.05 .15],[.15 .15]);

for regionName = REGIONS
    
    regionName = regionName{1};
    for dataInfo = [getDataInfo(regionName,pasPeriod),getDataInfo(regionName,futPeriod)]
        Imap = [];
        % try
        for ensNo = dataInfo.ensNo(1:end) % dataInfo.ensNo(:)'
            ensNo = ensNo{1};
            load(sprintf('%s%sIDAFs_%s_ensNo%s_%s.mat',dataInfo.savePath,filesep,...
                dataInfo.region.Name,ensNo,dataInfo.Years),'T','TIMES');
            pl = 0;% not plot
            parmhat_gev = grpstats(T.am,T.area,{@(x)fitAM(x,pl)});
            if size(parmhat_gev,1)==41
                parmhat_gev = [parmhat_gev(1:end,:);parmhat_gev(end,:)];
            end
            if ismatrix(parmhat_gev)
            parmhat_gev = mat2cell(parmhat_gev,ones(1,size(parmhat_gev,1)),[3]);
            end
            AREA = dataInfo.idaf.area';% grpstats(T.area,T.area,{'median'});
            RT = 1.1:0.1:40;
            I0 = getI(RT,parmhat_gev);
            [rt0,area0] = meshgrid(RT,AREA);
            Imap = cat(3,Imap,I0);
        end
        % catch me
        %     disp(ensNo);
        % end
        
        if strcmp(dataInfo.Years,pasPeriod)
            pasImap = Imap;
        elseif strcmp(dataInfo.Years,futPeriod)
            futImap = Imap;
        end
        
%         I = nanmedian(Imap,3);
%         I = griddata(rt0,area0,I,figForm.rt,figForm.area);
%         figure;setFigureProperty('Paper')
%         hold on;
%         pcolor(rt,log2(area),I);shading flat;
%         [h,c] = contour(rt,log2(area),I,[20:20:100],'w-','ShowText','on');
%         clabel(h,c,'Color',[0.9 0.9 0.9])
%         c.LineWidth = 1;
%         cptcmap('cw1-002','mapping','scaled','ncol',100);
%         set(gca,'YScale','linear')
%         ax = gca;
%         ax.YLim = log2([5,8000]);
%         ax.YTick = log2([5,10,20,50,100,200,500,1000,2000,4000,8000]);
%         ax.YTickLabels = 2.^ax.YTick;
%         ax.TickDir = 'out';
%         caxis([0,150]);
%         colorbar
%         title(sprintf('%s-%s',regionName,dataInfo.Years))
%         xlabel('Return Period (yr)')
%         ylabel('Catchment Area (km^2)')
%         savePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\IDAF\';
%         saveName = sprintf('%s-IDAF-%s-%s-Dur%02dh',regionName,dataInfo.Years,'allEns',1);
%         savePlot([savePath,filesep,saveName],'units','points','XYWH',[150,0,300,230],'needreply','N','onlyPng',true);
       
    end
    
    %%
    
    prcval = 50;
    Z = (futImap-pasImap)./pasImap;% Z = nanmedian(futImap-pasImap,3);
    Zs{itag} = Z;
    Z = 100*prctile(Z,prcval,3);%nanmedian(Z,3);
    
    sigZ = testSig(permute(futImap,[3,1,2]),permute(pasImap,[3,1,2]));
    Z = griddata(rt0,area0,Z,figForm.rt,figForm.area);
    sigZ = griddata(rt0,area0,sigZ,figForm.rt,figForm.area,'nearest');
    
    subplot(1,3,itag);
    % axes(ha(itag))
    setFigureProperty('Paper')
    set(0,'defaultAxesFontSize',10);
    % hold on;
    % Z = conv2(Z,ones(2,7)/14,'same');
    ff = @(x)imresize(x,'Scale',3,'method','bilinear');
    [rt0,area0] = meshgrid(rt,area);
    pcolor(ff(rt0),log2(ff(area0)),ff(Z));shading flat;hold on;
    % alpha(0.95)%'diff_darkBlue_darkRed'%BlueDarkRed18
    cptcmap('RdYlBu_11','mapping','scaled','ncol',8,'flip',true);% ncol 20
    % [~,c] = contour(rt,log2(area),Z,[-30,-20,-10,10,15,20,25,30],'k-','ShowText','on');
    [h,c] = contour(ff(rt0),log2(ff(area0)),ff(Z),[-40,-20,-10,10,20,30],'--',...
        'LineColor',[0.95 0.95 0.95],'ShowText','on');
    clabel(h,c,'Color',[0.9 0.9 0.9])
    [index] = find(Z > prctile(Z(:),99.5));
   %  plot(rt0(index),log2(area0(index)),'rx','markersize',10);hold on;
    plot(figForm.rt(mod(figForm.rt,2)>=0 & sigZ==0),...
        log2(figForm.area(mod(figForm.rt,2)>=0 & sigZ==0)),'.','color',[0.2 0.2 0.2 0.5],...
        'markersize',2);
    c.LineWidth = 1;
    ax = gca;
    ax.YScale = 'linear';
    ax.XScale = 'log';
    ax.YLim = log2([6,8000]);
    ax.XLim = [2,39.9];
    ax.YTick = log2([10,20,50,100,200,500,1000,3000,8000]);
    ax.YTickLabels = 2.^ax.YTick;
    ax.XTick = [2,5,10,20,40];
    ax.TickDir = 'out';
    % caxis([-35,35]);
    caxis([-40,40]);
    box off
    axis('square');
    % c = colorbar;
    % c.Ticks = [-50,-30,-10,0,10,30,50];
    % c.TickLabels = strcat(c.TickLabels,'%');
    % title(sprintf('(future-past)/past'))
    
    text(35,12.2,getPrintName(regionName),...
        'fontsize',10,'backgroundcolor',[0.9,0.9,0.9,0.1],...
        'horizontalalignment','right')
    
    xlabel('T (yr)')
    if itag == 1
    ylabel('A (km^2)')
    end
    savePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\IDAF\';
    
    itag = itag+1;
end

% format_xylabel(ha,1,3)

c = colorbar('EastOutside','position',[0.9 0.19 0.018 0.50]);
% c.Ticks = [-50,-30,-15,0,15,30,50];
c.Ticks = [-30,-20,-10,0,10,20,30];
c.TickLabels = strcat(c.TickLabels,'%');

% saveName = sprintf('%s-diff_IDAF-%s-Dur%02dh_percentile%02d',regionName,'allEns',1,prcval);
saveName = sprintf('diff_%s_IDAF-%s-Dur%02dh_percentile%02d','cpmRegions','allEns',1,prcval);
savePlot([savePath,filesep,saveName],'units','centimeters',...% 'XYWH',[5,0,18,7]
    'targetSize','dc','needreply','Y','onlyPng',false);
close all
%%
% easier-understanding plot
figure;
setFigureProperty('Paper');
set(0,'defaultAxesFontSize',10);
hold on;
h = handle(0);
ax = gca;% xlim([2,(8000)])
setXTick(ax);
col = {[66 146 199]/255,[240  60  43]/255,[0.3 0.3 0.3]};
for i = [1,3,2] %regionName = {'SCO','WAL','EUK'}
    subplot(1,3,i);hold on;
    Zp50 = 100*prctile(Zs{i},50,3);%
    Z0 = Zp50.*(repmat(200./RT,size(Zp50,1),1))/nansum(200./RT);
    h(i) = plot(log2(AREA),smooth(nansum(Z0,2),1),'-.','Color',col{i},'linewidth',1.5);
    
    % Zmean = 100*nanmean(Zs{i},3);
    % Z0_mean = Zmean.*(repmat(200./RT,size(Zmean,1),1))/nansum(200./RT);
    % plot(log2(AREA),nansum(Z0_mean,2),'o','Color',col{i},'linewidth',1);
    % close all

    Zp20 = 100*prctile(Zs{i},10,3);%
    Z10 = Zp20.*(repmat(200./RT,size(Zp20,1),1))/nansum(200./RT);
    Zp80 = 100*prctile(Zs{i},90,3);%
    Z90 = Zp80.*(repmat(200./RT,size(Zp80,1),1))/nansum(200./RT);
    fill([log2(AREA);flip(log2(AREA))],[nansum(Z90,2);flip(nansum(Z10,2))],h(i).Color,'LineStyle','none');
    alpha(0.12/i)
    axis('square');
    setXTick(gca);
%     for perval = linspace(1,100,12)
%     Zp_1 = 100*prctile(Zs{i},perval,3);%
%     Z_1 = Zp_1.*(repmat(200./RT,size(Zp_1,1),1))/nansum(200./RT);
%     plot(log2(AREA),nansum(Z_1,2),'Color',h(i).Color,'linewidth',0.5)
%     end
    % close all
    
    % SNR
    Z_12 = squeeze(nansum(Zs{i}.*(repmat(200./RT,...
        size(Zp20,1),1))/nansum(200./RT),2));
    SNR = nanmean(Z_12,2)./std(Z_12,[],2);
    SNR = table(SNR,AREA);
    
    
end
% title('Change of I(A,1h,RT)')

% box on;

c = legend(h,getPrintName(REGIONS),'Position',[0.03,0.2,0.4,0.04]);
c.NumColumns = 3;c.Box = 'off';

savePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\IDAF\';
saveName = sprintf('%s_ExtremeChange(CatchmentArea)','cpmRegions');
savePlot([savePath,filesep,saveName],'units','centimeters',...
    'targetSize','dc','needreply','Y','onlyPng',false);% 'XYWH',[5,0,18,7],

%% Auxillary function
function setXTick(ax)
ax.XLim = log2([5.3,8000]);
ax.YLim = [-10,51];
ax.XTick = log2([5,12,25,50,100,200,500,1000,2000,4000,8000]);
ax.XTickLabels = 2.^ax.XTick;
xlabel('A /Km^2');
ylabel('Relative Change(%)');
set(gca,'TickDir','Out')
xtickangle(90)
end
function [isSigDif] = testSig(Z1,Z2)
% Z: [loc1,loc2,sampleSize]

isSigDif = squeeze(ttest(Z1,Z2,'alpha',0.1));

end
function [pname] = getPrintName(REGIONS)
if ischar(REGIONS)
    switch(REGIONS)
        case 'CPM_NW'
            pname = 'NW UK';
        case 'CPM_NE'
            pname = 'NE UK';
        case 'CPM_S'
            pname = 'S UK';
    end
elseif iscell(REGIONS)
    pname = [];
    for i = 1:length(REGIONS)
        pname = cat(2,pname,{getPrintName(REGIONS{i})});
    end
end
end
function [I] = getI(RT,parmhat)
I = [];
P = 1-(1./RT);
for leni = 1:length(parmhat)
    try
    if ~isnan(parmhat{leni}(3))% GEV
        I(leni,:) = gevinv(P,parmhat{leni}(1),parmhat{leni}(2),parmhat{leni}(3));
    elseif  ~isnan(parmhat{leni}(1))% Gumbel
        I(leni,:) = evinv(P,parmhat{leni}(1),parmhat{leni}(2));
    else
        I(leni,:) = NaN(size(RT));
    end
    catch me
        1;
    end
end
end
function [parmhat] = fitAM(a,pl)
x = (.5:1:(numel(a)-.5))'./ numel(a);
r = 1./(1-x);
T = 1.05:0.05:40;
P = 1-(1./T);

% get stats for a series of a;
warning('')
options = statset('gevfit');
options.MaxIter = 500;
warning off
[parmhat,parmci] = gevfit(a,0.05,options);%%%%%%%%%%%%%%evfit(a);% 
X = gevinv(P,parmhat(1),parmhat(2),parmhat(3));%%%%%%%%%%%%%%%%%%%%%
[warnMsg,warnId] = lastwarn;
warning on
if ~isempty(warnMsg)
    warnMsg = '';
    [parmhat,parmci] = evfit(a);% 
    X = evinv(P,parmhat(1),parmhat(2));
    parmhat(3) = NaN;
    [warnMsg,warnId] = lastwarn;
end
warning off
if ~isempty(warnMsg)
parmhat(1:3) = NaN;
end

if pl
    colorI = 'r';
    hobs = plot(r,sort(a),'ko','markerfacecolor',colorI,'markersize',5); % cdf of y
    hold on;
    hfit = plot(T,X,'color',colorI);
    xlim([1,40])
    set(gca,'XScale','log')
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
idaf.area = [1:42].^2*(dataInfo.gridReso).^2;
dataInfo.idaf = idaf;

end

