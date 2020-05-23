function [LIB] = Main_OfflineTraining(version,eventNo4val,testConfig)
%
% This script will achieve OfflineTraining required for ANALOGUE
%
% programming paradigm: POP
%
% @ Yuting Chen
% 2020.02.24

arguments
    version (1,:) char = 'v1';
    eventNo4val (1,1) double = 1;
    testConfig (1,:) struct = struct();
end


tic

printInfo(version,eventNo4val);

libInfo = configureLibInfo(version,eventNo4val,testConfig);

RainLib = getRainfallLib(libInfo);

ClimateLib = getClimateLib(libInfo,RainLib);

FloodLib = getFloodLib(libInfo,RainLib.train.raintime,RainLib.test.raintime);

FCSLib = getFCSLib(libInfo);

LIB = struct('libInfo',libInfo,'RainLib',RainLib,'ClimateLib',ClimateLib,...
    'FloodLib',FloodLib,'FCSLib',FCSLib);

toc

end


%% AUXILLARY FUNC

function libInfo = configureLibInfo(version,eventNo4val,testConfig)

addpath(genpath('C:\Users\Yuting Chen\Dropbox (Personal)\MatLab_Func'));
addpath(genpath(cd))

libInfo = struct;

libInfo.climateDt = strcmp(version,'v1')*6+strcmp(version,'v2')*6;
libInfo.climatePriorTime = strcmp(version,'v1')*3+strcmp(version,'v2')*1;
libInfo.climatePressure = 700;% mustbemember{700,500}
libInfo.climateSmallerSpace = strcmp(version,'v1')*0+strcmp(version,'v2')*1;

% get time
libInfo.year = 2005:2017;

% get eventNo
if isfield(testConfig,'rainfallSource') && strcmp(testConfig.rainfallSource,'RAD')
    libInfo.rainfallSource = testConfig.rainfallSource;
    libInfo.rainfilePath = 'K:\DATA_FCS\RainEvent\Birm-34-30-Radar';
elseif (~isfield(testConfig,'rainfallSource')) || strcmp(testConfig.rainfallSource,'KED')
    libInfo.rainfallSource = testConfig.rainfallSource;
    libInfo.rainfilePath = 'K:\DATA_FCS\RainEvent\Birm-34-30';
end
FNS = dir([libInfo.rainfilePath,filesep,'EventNo*.mat']);
excluNo = 1:length(FNS);
libInfo.EventNoval = eventNo4val;
excluNo(eventNo4val) = [];
libInfo.EventNo = excluNo;

% get savePath
libInfo.savePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_FCS';

end

function RainLib = getRainfallLib(libInfo)
%
% GET ALL INFO which will be saved into RAINFALLLIB

RainLib = struct;
RainLib.train = struct;
RainLib.test = struct;
RainLib.EE = [];
RainLib.NN = [];

% IMPORT RAINFALL DATA: ALL NEEDED RADAR DATA

[RainLib.train.RIMAGES, RainLib.train.CRIMAGES, RainLib.train.raintime,...
    RainLib.test.RIMAGES, RainLib.test.CRIMAGES, RainLib.test.raintime,...
    RainLib.EE, RainLib.NN] = getRainEvents(...
    libInfo.year, libInfo.EventNo, libInfo.EventNoval,libInfo.rainfilePath);

[RainLib.train.PREDRAIN, RainLib.test.PREDRAIN] = formulateRainPredictor(...
    RainLib.train.RIMAGES, RainLib.test.RIMAGES,...
    RainLib.train.CRIMAGES, RainLib.test.CRIMAGES);


end

function ClimateLib = getClimateLib(libInfo,RainLib)

ClimateLib = struct;

pressure = libInfo.climatePressure;%700;
leadLagtime = libInfo.climatePriorTime;
DT = libInfo.climateDt;%desired reoslution of climate data

% rainfall event range: [390 419 270 303]
spaceEN = (libInfo.climateSmallerSpace==0)*zeros(1,4)+(libInfo.climateSmallerSpace==1)*[330 479 210 363];

% IMPORT ERA5 + Filter
[datainfo,E,N] = defineEra5DataBound(libInfo.year,pressure,spaceEN);
[CLIMATE,time,missval,dt] = getERA5(datainfo);

windowSize = DT/dt; % movmean(X,12*2,3)
b = (1/windowSize)*ones(1,windowSize);
CLIMATE0 = structfun(@(X)filter(b,1,X,[],3),CLIMATE, 'UniformOutput', false);% MA filter
% figure;
% plot(squeeze(CLIMATE0.T(1,1,:)),'r');hold on
% plot(squeeze(CLIMATE.T(1,1,:)),'k');hold on
% legend('Filtered Data','Input Data');
CLIMATE = CLIMATE0;


% IMPORT CORRESPONDING CLIMATE DATA: ALL NEEDED ERA5 DATA
[ClimateLib.train,ClimateLib.test] = formulateClimateData(...
    RainLib.train.raintime,...
    RainLib.test.raintime, RainLib.EE, RainLib.NN, ...
    CLIMATE, E, N, time, missval, leadLagtime);

ClimateLib.E = E;
ClimateLib.N = N;
ClimateLib.train.realtime = RainLib.train.raintime - leadLagtime;
ClimateLib.test.realtime = RainLib.test.raintime - leadLagtime;

end

function FloodLib = getFloodLib(libInfo,traintime,testtime)

FloodLib = struct;

[FloodLib.train, FloodLib.test] = deal(struct);

[hist_startTime, hist_endTime, hist_eventNos] = getHistEventTime(libInfo.rainfallSource);
[FloodLib.train.eventStats, FloodLib.test.eventStats,...
    FloodLib.train.eventNo, FloodLib.test.eventNo,...
    FloodLib.train.FilePath, FloodLib.test.FilePath] = deal([]);

for thistime = traintime'
    try
        FloodLib.train.eventNo(end+1) = findEventNo(thistime,hist_startTime,hist_endTime,hist_eventNos);
        FloodLib.train.FilePath{end+1} = ['G:\BIGDATA\TOPIC 2\ProcessedFiles\FloodMaps_40_30\',...
            sprintf('FloodMaps_Merged_EventNo%03d.mat',FloodLib.train.eventNo(end))];
        FloodLib.VarName = 'FMap';
    catch me
        me;
    end
end

for thistime = testtime'
    
    FloodLib.test.eventNo(end+1) = findEventNo(thistime,hist_startTime,hist_endTime,hist_eventNos);
    FloodLib.test.FilePath{end+1} = ['G:\BIGDATA\TOPIC 2\ProcessedFiles\FloodMaps_40_30\',...
        sprintf('FloodMaps_Merged_EventNo%03d.mat',FloodLib.test.eventNo(end))];
    FloodLib.VarName = 'FMap';
    
end


% get FloodStatistics
for i = 1:length(FloodLib.train.eventNo)
    
    FloodLib.train.eventStats(i) = getStats(FloodLib.train.FilePath{i},FloodLib.VarName);
    
end
for i = 1:length(FloodLib.test.eventNo)
    
    FloodLib.test.eventStats(i) = getStats(FloodLib.test.FilePath{i},FloodLib.VarName);
    
end


    function eventNo = findEventNo(time0,startTime,endTime,eventNos)
        eventNo = NaN(size(time0));
        for evi = 1:length(time0)
            eventNo(evi) = eventNos(startTime<=datenum(time0(evi)) & endTime>=datenum(time0(evi)));
        end
    end

    function stats = getStats(filename,filevar)
        FMAP = getfield(load(filename,filevar),filevar);
        stats = nansum(FMAP(:)>5);
        
    end
end

function FCSLib = getFCSLib(libInfo)

FCSLib = struct;

end




function printInfo(version,eventNo4val)
switch (version)
    case 'v1'
        fprintf('V1 training was used.\n');
        fprintf('Basic Information:\n')
        fprintf('Cliamte Aggregation duration: 6 hour\n')
        fprintf('EventNo %3d excluded.\n',eventNo4val)
    case 'v2'
        fprintf('V2 training was used.\n');
        fprintf('Basic Information:\n')
        fprintf('Cliamte Aggregation duration: 6 hour\n')
        fprintf('Climate Area Only Birm included (1h before).\n')
        fprintf('EventNo %3d excluded.\n',eventNo4val)
    otherwise
        %
        
end


end






%% AUX 2
function [PREDRAIN,PREDRAINval] = formulateRainPredictor(RIMAGES,RIMAGESval,CRIMAGES,CRIMAGESval)

PREDRAIN = struct;
RVec = reshape(RIMAGES,numel(RIMAGES(:,:,1)),[]);
PREDRAIN.WAR = nansum(RVec>0.01,1)./numel(RVec(:,1));
PREDRAIN.IMF = nanmean(RVec,1);%./PREDRAIN.WAR;
PREDRAIN.MAX = nanmax(RVec,[],1);


PREDRAINval = struct;
RVec = reshape(RIMAGESval,numel(RIMAGESval(:,:,1)),[]);
PREDRAINval.WAR = nansum(RVec>0.01,1)./numel(RVec(:,1));
PREDRAINval.IMF = nanmean(RVec,1);%./PREDRAINval.WAR;
PREDRAINval.MAX = nanmax(RVec,[],1);

end


function [TRAIN,TEST] = formulateClimateData(raintime,raintimeval,EE,NN,CLIMATE,E,N,time,missval,leadLagtime)

% CLIMATE: is <struct> including: ('GEOZ',GEOZ,'O3',O3,'RH',RH,'CRWC',CRWC,'T',T,'U',U,'V',V);

[ni,nlen] = findInd(NN(:,1),N);
[ei,elen] = findInd(EE(1,:),E);

% leadLagtime;% need to change. because it refers to the distance to peak time not the current time.
% LEDLAGtIME:3; % ref: (Loriaux et al., 2016) P5481

% locate only birm area, also to look at the into <$leadLagtime$> hours before
tind = dsearchn(datenum(time),datenum(raintime))-leadLagtime;
TRAIN = structfun(@(X)X(ni:ni+nlen-1,ei:ei+elen-1,tind),CLIMATE, 'UniformOutput', false);

tind = dsearchn(datenum(time),datenum(raintimeval))-leadLagtime;
TEST = structfun(@(X)X(ni:ni+nlen-1,ei:ei+elen-1,tind),CLIMATE, 'UniformOutput', false);

% reshape and get only one val for each snapshot
TRAIN = structfun(@(X)reshape(X,nlen*elen,[]),TRAIN, 'UniformOutput', false);
TEST = structfun(@(X)reshape(X,nlen*elen,[]),TEST, 'UniformOutput', false);

TRAIN = structfun(@(X)getMean(X,missval,1),TRAIN, 'UniformOutput', false);
TEST = structfun(@(X)getMean(X,missval,1),TEST, 'UniformOutput', false);

    function XM = getMean(X,missval,meanDim)
        X(X==missval) = NaN;
        XM = nanmean(X,meanDim);
    end

    function [ni,nlen] = findInd(Vec,Vwhole)
        nimax = dsearchn(Vwhole,max(Vec)+14);
        nimin = dsearchn(Vwhole,min(Vec)-14);
        nlen = abs(nimin-nimax)+1;
        ni = min(nimin,nimax);
    end

end


function [datainfo,E,N] = defineEra5DataBound(year,pressure,spaceEN)

datainfo = struct;

filePath = 'K:\DATA_FCS\ERA5_Birm\ERA5_pressure';
fileName = sprintf('era5_%dhPa_Birmingham_all_1979_2018.nc',pressure);

datainfo.fileN = [filePath,filesep,fileName];
A = ncinfo(datainfo.fileN);

% ADDITIONAL GRID
% fileName = sprintf('era5_%dhPa_Birmingham_all_additional_1979_2018.nc',pressure);
% datainfo.fileN = [filePath,filesep,fileName];
% Aadd = ncinfo(datainfo.fileN = [filePath,filesep,fileName];);

% For near surface variables
filePath = 'K:\DATA_FCS\ERA5_Birm\ERA5_singleLevel';
fileName = sprintf('era5_singleLevel_all_Birmingham_1979_2018.nc');
datainfo.fileN_singleLevel = [filePath,filesep,fileName];
A = ncinfo(datainfo.fileN_singleLevel);

fprintf(['Notice that for full datasets around $Birmingham$',...
    'might need to include $Add$\n'])

LON = ncread(datainfo.fileN,'longitude');
LAT = ncread(datainfo.fileN,'latitude');
[E, N] = ll2os(LAT, LON);
E = E/1000;
N = N/1000;
time = ncread(datainfo.fileN,'time')/24+datenum(datetime(1900,1,1,0,0,0));
datainfo.timestart = find(datetime(datevec(time)).Year==year(1),1);
datainfo.timeend = find(datetime(datevec(time)).Year==year(end)+1,1)-1;
datainfo.t_len = datainfo.timeend-datainfo.timestart+1;

if ~any(spaceEN==0)
    datainfo.loni = find(min(abs(E-spaceEN(1))) == abs(E-spaceEN(1)));
    datainfo.lonlen = find(min(abs(E-spaceEN(2))) == abs(E-spaceEN(2)))-datainfo.loni+1;
    datainfo.lati = find(min(abs(N-spaceEN(4))) == abs(N-spaceEN(4)));
    datainfo.latlen = find(min(abs(N-spaceEN(3))) == abs(N-spaceEN(3)))-datainfo.lati+1;
else
    datainfo.loni = 1;
    datainfo.lonlen = Inf;
    datainfo.lati = 1;
    datainfo.latlen = Inf;
end


% Result Checking Module
if any([datainfo.loni,datainfo.lonlen,datainfo.lati,datainfo.latlen]<=0)
    f = msgbox('Climate Data Format (N Axis) is changed! Please Contact Yuting to modify code', 'Error','error');
end

if any([E(:),N(:)] > 1e4)
    f = msgbox('Climate Data Format (E,N) is changed! Please Contact Yuting to modify code', 'Error','error');
end

end



