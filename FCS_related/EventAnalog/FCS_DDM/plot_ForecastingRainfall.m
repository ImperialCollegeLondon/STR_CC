% --------------------------------------------------------------- %
% THIS IS TO CHECK & PLOT ENSEMBLES_BASED FORECASTING MODEL
% --------------------------------------------------------------- %

leadT = 2; % unit [h]
ensNo = 12;
eventNo = 61;
expNo = 'C1P5-Rad1';
[~,~,statFp] = getExpNoInfo(expNo);

saveFigFp = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_FCS\ForecastingResTest\';
Area = getBoundaryShp('Birm');
mkdir(saveFigFp);

for evi = eventNo(:)
    
    statFn = [statFp,sprintf('STATS_v1_ev%03d.mat',evi)];
    [obsfp,obsid,simfp,simid,forid,flmaps] = getDat(statFn);
    Obs = load(obsfp,'EE','NN','RMap');
    
    for evti = 145 % round(size(simfp,1)*0.1):round(size(simfp,1)*0.8)
        
        figName = sprintf('Forecast%03d-%04d',evi,evti);
        close all
        
        plotAllinOne(Area,Obs.EE,Obs.NN,squeeze(Obs.RMap(:,:,evti)),...
            simfp(evti,:),simid(evti,:),forid(evti,:),...
            squeeze(flmaps.val(evti,:,:,:)),flmaps.E,flmaps.N);
        savePlot([saveFigFp,figName],'Units','centimeters','XYWH',[0,0,14,7.5],'needreply','N','onlyPng',true);
    
    end
    
end


%% AUX FUNC
function [obsfp,obsid,simfp,simid,forid,flmaps] = getDat(statFn)
% Output:
% obsfp char
% obsid (:,1) double
% simfp {:,1} <cell>.<char>
% simid (:,12) double

load(statFn,'FOREOUTPUT','STATS')
obsevi = FOREOUTPUT.obsEvNo(1);
obsfp = ['K:\DATA_FCS\RainEvent\Birm-34-30-Radar\',sprintf('EventNo%03d.mat',obsevi)];
obsid = STATS.eventTi+1;
simfp = arrayfun(@(x)['K:\DATA_FCS\RainEvent\Birm-34-30-Radar\',sprintf('EventNo%03d.mat',x)], ...
    FOREOUTPUT.analogPredEvNo, 'UniformOutput', false);
EventTable = readtable('K:\DATA_FCS\RainEvent\Birmingham_SelectedEvents_DataDriven_RadarVersion.csv');
simid = cell2mat(arrayfun(@(x)getEviTimeId(x,EventTable), FOREOUTPUT.analogTime, 'UniformOutput', false));
forid = simid + 12;

flmaps = struct('val',FOREOUTPUT.FLO_pred.floodmaps,'E',FOREOUTPUT.E/1000,'N',FOREOUTPUT.N/1000);

end


function [ind] = getEviTimeId(imageTime,T)

% imageTime = double(string(imageTime));
T.startTime = datetime(string(T.startTime),'Format','yyyyMMddHHmm');
T.endTime = datetime(string(T.endTime),'Format','yyyyMMddHHmm');
ind = round(minutes(imageTime-T.startTime(imageTime>=T.startTime & imageTime<=T.endTime))./5)+1;

if length(ind) ~=1 || ind <=0
    uiwait(msgbox('check time id','Forecasting','Searching'));
end

end

function plotAllinOne(Area,EE,NN,obsMap,simfp,simid,forid,fmap,fE,fN)

ha = tight_subplot(6,12,[0.05 0.02],[.05 .05],[.05 .05]);
view(2)
ha2 = tight_subplot(6,12,[0.05 0.02],[.05 .05],[.05 .05]);
set(ha([1:5,7:12]),'visible','off')
set(ha2(:),'visible','off')
axes(ha(6))
plotThis(EE,NN,obsMap)


for ensNo = 1:12
    load(simfp{ensNo},'RMap')
    axes(ha(12+ensNo))
    plotThis(EE,NN,squeeze(RMap(:,:,simid(ensNo))));
    
    axes(ha(24+ensNo))
    try
        plotThis(EE,NN,squeeze(RMap(:,:,forid(ensNo))));
    catch
        plotThis(EE,NN,EE*0);
    end
%     
%     axes(ha(36+ensNo))
%     try
%         plotThis(EE,NN,squeeze(RMap(:,:,forid(ensNo)+1)));
%     catch
%         plotThis(EE,NN,EE*0);
%     end
%     
%     axes(ha(48+ensNo))
%     try
%         plotThis(EE,NN,squeeze(RMap(:,:,forid(ensNo)+2)));
%     catch
%         plotThis(EE,NN,EE*0);
%     end
%     
%     axes(ha(60+ensNo))
%     try
%         plotThis(EE,NN,squeeze(RMap(:,:,forid(ensNo)+3)));
%     catch
%         plotThis(EE,NN,EE*0);
%     end
    
    
    
    for jj = 1
    if 1
    ax2 = ha2(36+ensNo);
    ax1 = ha(36+ensNo);
    axes(ha(36+ensNo))
    plotBirmBackground100(ha(36+ensNo))
    hold on;
    axes(ha2(36+ensNo))
    thismap = squeeze(fmap(:,:,ensNo));
    thismap(thismap == 0) = NaN;
    [cmap, ~, ~, ~, ~] = cptcmap('flood_blue','mapping','scaled','ncol',15);
    colfmap = NaN(numel(thismap),3);
    colfmap(~isnan(thismap(:)),:) = cmap(getLevel(thismap(~isnan(thismap)),linspace(0,20,15)),:);
    scatter(ax2,fE(:),fN(:),thismap(:)*0.4,colfmap,'fill');shading flat;
    alpha(0.6)
    %%Link them together
    linkaxes([ax1,ax2])
    %%Hide the top axes
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
%     plotThis(fE,fN,);
%     cptcmap('flood_blue','mapping','scaled');
%     caxis([0,20])
    xlim([min(EE(:)),max(EE(:))])
    ylim([min(NN(:)),max(NN(:))])
    set(gca,'YTick',[],'YTickLabel',[],'XTick',[],'XTickLabel',[],'Linewidth',1);
    end
    end
end
    function plotThis(EE,NN,Rmap)
        Rmap(Rmap == 0) = NaN;
        pcolor(EE,NN,Rmap);
        shading flat;
        hold on;
        plot(Area.X/1000,Area.Y/1000,'k','linewidth',1);
        cptcmap('mld_rain-mmh', 'mapping','scaled');%,'ncol',20);
        caxis([0,80])
        % cptcmap('precip_meteoswiss', 'mapping','direct');
        % caxis([0,3])
        % axis off% precip_meteoswiss% cw1-002
        set(gca,'YTick',[],'YTickLabel',[],'XTick',[],'XTickLabel',[],'Linewidth',1);
        set(gca,'color',[0.8 0.8 0.8]);
    end
end


function plotBirmBackground100(ax1)
% ax1 = axes;
load('G:\BIGDATA\TOPIC 2\BirmLocalMap_reso100m.mat','x','y','imageData');

Area = getBoundaryShp('Birm');
XLIM = [Inf,0];YLIM = [Inf,0];
Area.X = Area.X/1000;
Area.Y = Area.Y/1000;
for i = 1:length(x)
    for tileNum = 1:numel(x{i})
        pcolor(x{i}{tileNum}/1000,y{i}{tileNum}/1000,imageData{i}{tileNum}); shading flat; hold on;
        XLIM = [min(XLIM(1),min(x{i}{tileNum})/1000),max(XLIM(2),max(x{i}{tileNum}/1000))];
        YLIM = [min(YLIM(1),min(y{i}{tileNum})/1000),max(YLIM(2),max(y{i}{tileNum})/1000)]; 
    end
    cptcmap('GMT_gray','mapping', 'scaled','flip',false,'ncol',256); caxis([0,255])
    drawnow;
end
plot(Area.X,Area.Y,'k','linewidth',1)
set(gca,'YDir','normal')
% axis equal
axis off
xlim(XLIM);ylim(YLIM)
clear x y imageData

end

function mapLevel = getLevel(maps,lowThre)
mapLevel = NaN(size(maps));
lowThre = [lowThre,inf];
for li = 1:numel(lowThre)-1
    mapLevel(maps>=lowThre(li) & maps<lowThre(li+1)) = li;
end
end

