% THIS SCRIPT IS TO EXTRACT THE RADAR DATA FOR THE PERIOD OF 2005-2017
%
% Following functions will be achieved.
% (1) Extracting whole dataset
% (2) Seperate event based on:
% (3) Possible merging??
%
% in order to provide input for classification scheme
%
% @ Yuting Chen
% Update: 2020/02/26
%
% Note:


clear;clc
close all

savePath = 'K:\DATA_FCS';


%% EXTRACTING WEATHER RADAR RAIN

load('K:\DATA_FCS\RainEvent\Birm-34-30\EventNo004.mat','EE','NN')

XX = EE*1000;
YY = NN*1000;

for YEAR = 0%2005:2017
    
    [PRS,status] = importNIMROD_P(XX,YY,YEAR);

    save([savePath,filesep,'RadarRain\Nimrod5min_year',sprintf('%04d.mat',YEAR)],'PRS','XX','YY');

end

%% CORRESPONDING RADAR IMAGES for each identified flood-inducing rainfall event

[floodStartT,floodEndT] = getAllFloodEvents();

savePath = 'K:\DATA_FCS\RainEvent\Birm-34-30-Radar';
mkdir(savePath)
kedPath = 'K:\DATA_FCS\RainEvent\Birm-34-30';
filenamepre = [savePath,filesep,'EventNo'];

resoMin = 5;
YEAR = 0;
softT = 3;% unit:[hour]

[PRS,T] = deal([]);
% extract more than current extract event time using RGs.
for YEAR = 2005:2017
    T0 = datetime(YEAR,1,1,'format','default'):minutes(resoMin):...
        datetime(YEAR+1,1,1,'format','default')-minutes(resoMin);
    A = load(['K:\DATA_FCS\RadarRain\Nimrod5min_year',sprintf('%04d.mat',YEAR)],'PRS');
    A.PRS(isnan(A.PRS)|(A.PRS<0)) = 0;
    % A.PRS = double(A.PRS)/32;
    PRS = cat(3,PRS,A.PRS);
    T = [T,T0];
end
size(PRS)
size(T)
for evi = 1:length(floodStartT)
    floodStartT(evi) = floodStartT(evi)-hours(softT);
    floodEndT(evi) = floodEndT(evi)+hours(softT);
    tt_temp = floodStartT(evi):minutes(resoMin):floodEndT(evi);
    ti = find(abs(datenum(T-floodStartT(evi)))<1e-5):find(abs(datenum(T-floodEndT(evi)))<1e-5);
    if length(tt_temp)~=length(ti)
        1;
    end
    RMap = double(PRS(:,:,ti))/32;
    if isempty(RMap)
        1;
    end
    filename = [filenamepre,sprintf('%03d.mat',evi)];
    load(['K:\DATA_FCS\RainEvent\Birm-34-30\EventNo',sprintf('%03d.mat',evi)],'EE','EventTimes','NN')
    save(filename,'EE','EventTimes','NN','RMap');

end
clear PRS A

% EXCLUDE THOSE STARTING/ENDING PERIOD having low imf(wet grids) (<0.03) & war (<0.02)
savePath = 'K:\DATA_FCS\RainEvent\Birm-34-30-Radar';
filenamepre = [savePath,filesep,'EventNo'];
set(gca,'linewidth',2)

[startTime,endTime] = deal(string([]));
eventNos = [];
for evi = 1:length(floodStartT)
    
    filename = [filenamepre,sprintf('%03d.mat',evi)];
    load(filename,'EE','EventTimes','NN','RMap');
    R = reshape(RMap,prod(size(RMap,[1,2])),[]);
    WAR = nanmean((R>0.03),1);
    bar(WAR);
    title(sprintf('Event%03d',evi));
    ylabel('WAR')
    
    % for begining part.
    barWAR = cumsum(WAR)./(1:length(WAR));
    i = find(barWAR<0.02 & WAR<0.02);
    if ~isempty(i)
        i = i(end);
        hold on;plot(i,0,'ro','markersize',10);hold off;
    else
        i = 0;
    end
    
    % for ending part.
    barWAR = cumsum(flip(WAR))./(1:length(WAR));
    j = find(barWAR<0.02 & flip(WAR)<0.02);
    if ~isempty(j)
        j = j(end);
        hold on;plot(length(WAR)-j+1,0,'bo','markersize',10);hold off;
    else
        j = 0;
    end
   
    drawnow; % pause(0.5)

    timest = datetime(floodStartT(evi)+minutes(i*5),'format','yyyyMMddHHmm');
    
    RMap = RMap(:,:,i+1:length(WAR)-j);
    timese = datetime(timest+minutes((size(RMap,3)-1)*resoMin),'format','yyyyMMddHHmm');
    EventTimes = [char(timest),'-',char(timese)];

    save(filename,'EE','EventTimes','NN','RMap');
    startTime(evi,1) = string(timest);
    endTime(evi,1) = string(timese);
    eventNos(evi,1) = evi;
end

A = table(eventNos,startTime,endTime);
writetable(A,'K:\DATA_FCS\RainEvent\Birmingham_SelectedEvents_DataDriven_RadarVersion.csv');





%% EVENT SEPARATION

%{
tagEventNo = 1;

for YEAR = 2005:2017
    
    savePath = 'K:\DATA_FCS';
    load([savePath,filesep,'RadarRain\Nimrod5min_year',sprintf('%04d.mat',YEAR)],'PRS','XX','YY');
    PRS(isnan(PRS)) = 0;
    PRS = double(PRS)/32;
    
    %% SEPERATE EVENT
    
    % those events used in flood prediction will be taken out first.
    % and then other events will be identified and extractd
    
    % Seperation criteria:
    % which is the criteria used in
    % https://doi.org/10.1002/2013WR014437
    % Including:
    % compute WAR (a threshold of 0.1mm/h was used)
    % 1h MA window
    % at least 2h war <%2
    
    resoMin = 5;
    
    R = reshape(PRS,prod(size(PRS,[1,2])),[]);
    T = datetime(YEAR,1,1,'format','default'):1/24/(60/resoMin):...
        datetime(YEAR+1,1,1,'format','default')-1/24/(60/resoMin);
    WAR = nanmean((R>0.03),1);% a threshold of 0.03mm/h was used to compute the WAR
    windowSize = 60/resoMin;
    b = (1/windowSize)*ones(1,windowSize);
    WAR = filter(b,1,WAR);
    % WAR = smooth(WAR,(60/resoMin));% 1h MA window
    plot(WAR)
    figure;
    cdfplot(WAR)
    WAR0 = WAR;
    WAR(WAR <= 0.02) = 0;
    [lengthi,stormsi] = seperateRadarEvents( WAR, 0.02, 6, 60);
    
    
    close all
    
    figure;
    hold on
    bar((1:length(WAR))/12,WAR0);
    plot(stormsi/12,(stormsi>0)-1,'r.','markersize',50);
    plot((stormsi+lengthi-1)/12,(stormsi>0)-1,'k.','markersize',50);
    xlabel('time(hr)');ylabel('WAR(-)');
    set(gca,'linewidth',2);box on
    ax = gca;
    plot([ax.XLim],[0.02,0.02],'k--');
    ax.YLim = [0,1];
    fprintf('max duration(hr) in year %4d: %3d hr\n',YEAR,round(max(lengthi)/12))
    
    % figPath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_FCS';
    % filename = [figPath,filesep,'StormEventSeparation'];
    % savePlot(filename,'XYWH',[150,0,800,300],'needreply','Y');
    
    close all
    
    %% SAVE EACH EVENT
    
    stormst = T(stormsi);stormst = stormst';
    savePath = 'K:\DATA_FCS\RadarEvent';
    filenamepre = [savePath,filesep,'EventNo'];
    
    for evi = 1:length(stormsi)
        
        filename = [filenamepre,sprintf('%04d.mat',tagEventNo)];
        trange = stormsi(evi):(stormsi(evi)+lengthi(evi)-1);
        [EE,NN,EventTimes,RMap] = createEventInfo(stormst(evi),lengthi(evi),PRS(:,:,trange),XX,YY);
        save(filename,'EE','EventTimes','NN','RMap');
        tagEventNo = tagEventNo+1;
        
    end
    
end

%}
% POSSIBLE MERGING


%% GET SEVERAL FEATURE MAP
filepath = 'K:\DATA_FCS\RadarEvent\EventNo*.mat';
FNs = dir(filepath);


[floodStartT,floodEndT] = getAllFloodEvents();
[PRD5min,PRD1h,ARD,D1h] = deal(NaN(34,30,0));
FLOOD = [];
for fn = FNs'
    
    load([fn.folder,filesep,fn.name],'RMap','EventTimes');
    CRain_aux = squeeze(nansum(RMap,3)/12);%[mm]
    RMap1h = imresize3(RMap,'scale',[1,1,1/12],'method','box');
    PRain5min_aux = squeeze(nanmax(RMap,[],3));%[mm/h]
    PRain1h_aux = squeeze(nanmax(RMap1h,[],3));%[mm/h]
    floodId_aux = classifyFlood(EventTimes,floodStartT,floodEndT);
    
    PRD5min = cat(3,PRD5min,PRain5min_aux);
    PRD1h = cat(3,PRD1h,PRain1h_aux);
    ARD = cat(3,ARD,CRain_aux);
    FLOOD = [FLOOD,floodId_aux];
    
end

%% Check indicator:
figure;
boxplot(squeeze(nanmean(nanmean(ARD,1),2)),FLOOD,'whisker',100);



%%
% ANALYSING 2 GET CRITICAL DEPTH

filepath = 'K:\DATA_FCS\RadarEvent\EventNo*.mat';
FNs = dir(filepath);


[floodStartT,floodEndT] = getAllFloodEvents();% KED version

[x,x0,y,y0] = deal([]);

for fn = FNs'
    
    load([fn.folder,filesep,fn.name],'RMap','EventTimes');
    CRain = squeeze(nansum(RMap,3)/12);
    CRain = nanmax(CRain(:));
    if classifyFlood(EventTimes,floodStartT,floodEndT) == false
        x0 = [x0;round(size(RMap,3)/12)];
        y0 = [y0;round(CRain)];
    else
        x = [x;round(size(RMap,3)/12)];
        y = [y;round(CRain)];
    end
    
end

%%
figure;
setFigureProperty('Subplot2');
scatter(x0,y0,100,'k','filled');
hold on;
plot(x,y,'ro','markersize',5,'markerfacecolor','r');alpha(0.1)
box on;
set(gca,'linewidth',2)
ylim([0,150])
xlim([0,150])
legend('NoFlood Event 2005-2010','Flood Event 2005-2010');

xlabel('Duration (hour)');
ylabel('Cummulative Rainfall Depth (mm)');

figPath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_FCS';
filename = [figPath,filesep,'StormDD_Scatter'];
savePlot(filename,'XYWH',[150,0,800,300],'needreply','Y','onlyPng',true);

%%
figure;
setFigureProperty('Subplot2');
boxplot(y0,x0,'whisker',200);
hold on;
plot(x,y,'ro','markersize',5,'markerfacecolor','r');alpha(0.5)
box on;
ylim([0,80])
xlim([0,40])
set(gca,'linewidth',2)
xtickangle(90)
xlabel('Duration (hour)');
ylabel('Cummulative Rainfall Depth (mm)');
figPath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_FCS';
filename = [figPath,filesep,'StormDD_Boxplot'];
savePlot(filename,'XYWH',[150,0,800,300],'needreply','Y','onlyPng',true);

%%
% plot flooded

filepath = 'K:\DATA_FCS\RainEvent\Birm-34-30\EventNo*.mat';
FNs = dir(filepath);
[x,y] = deal([]);
for fn = FNs'
    load([fn.folder,filesep,fn.name],'RMap');
    CRain = squeeze(nansum(RMap,3)/12);
    CRain = nanmean(CRain(:));
    x = [x;round(size(RMap,3)/12)];
    y = [y;round(CRain)];
    
    
end
hold on;scatter(x,y,50,'g','filled');alpha(0.5)


%% AUXILLARY FUNCTION

function [floodStartT,floodEndT] = getAllFloodEvents()

[startTime,endTime,eventNos] = getHistEventTime('KED');

floodStartT = datetime(datevec(startTime));
floodEndT = datetime(datevec(endTime));

end




function binRes = classifyFlood(timeStr,floodStartT,floodEndT)


timeVec = getTimeVec(timeStr);
timeVec1 = datenum(timeVec(1));
timeVec2 = datenum(timeVec(end));

[ii] = find(abs(timeVec1-datenum(floodStartT))<12/24 & (timeVec2-datenum(floodEndT))<12/24);

binRes = ~(isempty(ii));

    function timeVec = getTimeVec(timeStr)
        % '201407081015-201407082000'
        timeVec = datetime(timeStr(1:12),'format','yyyyMMddHHmm'):1/24/12:...
            datetime(timeStr(14:end),'format','yyyyMMddHHmm');
        timeVec = reshape(timeVec,[],1);

    end

end


function [EE,NN,EventTimes,RMap] = createEventInfo(stormst,lengthi,PRS,XX,YY)

timest = char(datetime(stormst,'format','yyyyMMddHHmm'));
timese = char(datetime(stormst+lengthi/24/12,'format','yyyyMMddHHmm'));
EventTimes = [timest,'-',timese];

EE = round(XX/1000);
NN = round(YY/1000);

RMap = (NN(2)-NN(1)>0)*PRS + (NN(2)-NN(1)<=0)*flip(PRS,1);

end





