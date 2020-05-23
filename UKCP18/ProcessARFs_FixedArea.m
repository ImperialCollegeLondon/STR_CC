%
% This script is to calculate ARFs using a fixed-area approach.
%
% Data: 1980-2000 UKCP CPM rainfall or 2060-2080 UKCP CPM rainfall
%
%
%
% ref: Svensson, C.; Jones, D.A.. 2010 Review of methods for deriving areal
%     reduction factors. Journal of Flood Risk Management, 3. 232-245.
%     10.1111/j.1753-318X.2010.01075.
%

clear;clc
dataInfo = getDataInfo('WAL','1980-2000');

% % #for radar data#
% for ensNo = {'RAD'}
%     ensNo = ensNo{1};
%     for dataInfo = [getDataInfo('SCO','2007-2018'),...
%         getDataInfo('EUK','2007-2018'),...
%         getDataInfo('WAL','2007-2018')]

% % #for gear data#
% for ensNo = {'GEAR'}
%     ensNo = ensNo{1};
%     for dataInfo = [getDataInfo('SCO','1990-2014'),...
%         getDataInfo('EUK','1990-2014'),...
%         getDataInfo('WAL','1990-2014')]
    
% % #for cpm data#
for ensNo = dataInfo.ensNo(:)'
ensNo = ensNo{1};
for dataInfo = [getDataInfo('EUK','1980-2000'),getDataInfo('EUK','2060-2080'),...
        getDataInfo('SCO','1980-2000'),getDataInfo('SCO','2060-2080'),...
        getDataInfo('WAL','1980-2000'),getDataInfo('WAL','2060-2080')]

        [T,ALLRMAX,TIMES] = deal([]);
        
        % # cpm & radar        
        [data,time] = getData(dataInfo,ensNo,[]);
        for year = reshape(unique(time.Year),1,[])
            stats = []; % <table>
            thisYearIndex = time.Year == year;
            thisYearTime = time(thisYearIndex);
            thisYearData = data(:,:,thisYearIndex);
        % # cpm & radar     

%         % # gear1hr (is a bad dataset!)
%         for year = reshape(getYearRange(dataInfo),1,[])
%         stats = []; % <table>   
%         [thisYearData,thisYearTime] = getData(dataInfo,ensNo,year);    
%         % # gear
        
            for duri = 1:length(dataInfo.arf.duration)
                duration = dataInfo.arf.duration(duri);
                
                startT = tic;
                
                threshold = getThreshold(duration,dataInfo);
                data0 = processData(thisYearData,duration);% process data
                [allR_Max,allT_Max] = deal([]);
                cellRange = getCellRange(data0,dataInfo);
                
                [I_aMax] = getAreaMax(data0,dataInfo);
                for celli = cellRange
                    for cellj = cellRange
                        if ~ isscalar(celli)
                            error('check for loop for celli & cellj');
                        end
                        pData = getCurrentPD(data0,celli,cellj,dataInfo);
                        [R_pMax,R_aMax,T_pMax,T_aMax] = getMaxima(I_aMax,pData,thisYearTime,threshold,dataInfo);
                        [R_aMax,T_aMax] = onlyInArea(celli,cellj,size(data0),R_aMax,T_aMax,dataInfo);
                        % save this snap to allSnaps
                        if ~isempty(R_pMax)
                            allR_Max{length(allR_Max)+1,1} = [reshape(R_pMax,1,[]);reshape(R_aMax,1,[])];
                            allT_Max{length(allT_Max)+1,1} = [reshape(T_pMax,1,[]);reshape(T_aMax,1,[])];
                        end
                    end
                end
                %
                allR_Max = cellfun(@(x)reshape(x,[1,size(x)]),allR_Max,'UniformOutput',false);
                allR_Max = cell2mat(allR_Max);
                allT_Max = cat(1,[],allT_Max);
                
                [allR_Max,allT_Max] = filterOutSameStorm(allR_Max,allT_Max);
                
                am = allR_Max(:,1);
                arf = getArf(allR_Max,allT_Max,dataInfo);
                dur = repmat(duration,size(allT_Max));
                ens = repmat(ensNo,size(allT_Max));
                T = [T;table(allT_Max,dur,arf,ens,am)];
                ALLRMAX = cat(2,ALLRMAX,allR_Max(:)');
                TIMES = cat(2,TIMES,allT_Max(:)');
                fprintf('ensNo %s duration %02d %03.1f seconds \n',ensNo,duration,toc(startT));
            end
            
        end
        ALLRMAX = ALLRMAX';TIMES = TIMES';
        save(sprintf('%s%sARFs_%s_ensNo%s_%s.mat',dataInfo.savePath,filesep,...
            dataInfo.region.Name,ensNo,dataInfo.Years),'T','TIMES','ALLRMAX','-v7.3');
    end
    
end

%% AUXILLARY FUNCTION
function thisarf = getArf(allR_Max,allT_Max,dataInfo) %#ok<INUSL>
thisarf = NaN(size(allR_Max,1),length(dataInfo.arf.area));
for ii = 1:size(allR_Max,1)
    thisMax = squeeze(allR_Max(ii,:,:));
    tag = 0;
    for area = reshape(dataInfo.arf.area,1,[])
        tag = tag+1;
        thisarf(ii,tag) = getThisARF(thisMax(:,tag),area); 
    end
end

    function thisval = getThisARF(thisMax,area) %#ok<INUSD>
        thisval = thisMax(2)./thisMax(1);
    end
   
end
function threshold = getThreshold(duration,dataInfo)
threshold = dataInfo.arf.threshold(dataInfo.arf.duration == duration);
end

function [allR_Max,allT_Max] = filterOutSameStorm(allR_Max,allT_Max)
% no filter for Fixed-area method.
end
function [R_aMax,T_aMax] = onlyInArea(celli,cellj,size_data0,R_aMax,T_aMax,dataInfo)
tag = 0;
for area = reshape(dataInfo.arf.area,1,[])
    tag = tag+1;
    centi = (1+size_data0(1))/2;
    if (celli-(centi))^2+(cellj-centi)^2>(area/pi)
        R_aMax(tag) = NaN;
        T_aMax(tag) = datetime(0,0,0); 
    end
end

end
function [I_aMax] = getAreaMax(data_ori,dataInfo)
% return index
[I_aMax] = deal([]);
data0 = imresize3(data_ori,'Scale',[2.2,2.2,1],'method','box');
for area = 2.2*2.2*reshape(dataInfo.arf.area,1,[])
    % circle area
    [ii,jj] = meshgrid(1:size(data0,1),1:size(data0,1));
    centi = (1+size(data0,1))/2;
    [notCir] = ((ii-(centi)).^2+(jj-centi).^2)>(area/pi);
    curData_temp = reshape(data0,prod(size(data0,[1,2])),[]);
    curData_temp(notCir,:) = NaN;
    centralTS = squeeze(nansum(curData_temp,1));
    % plot(centralTS);
    [~,I] = sort(centralTS,'descend');
    I_aMax = [I_aMax,I(1)];
end

end
function [R_pMax,R_aMax,T_pMax,T_aMax] = getMaxima(I_aMax,pData,thisYearTime,threshold,dataInfo)

thisImageIndex = 1;

centralTS = pData;
[~,I] = sort(centralTS,'descend');
R_pMax = centralTS(I(thisImageIndex));
T_pMax = thisYearTime(I(thisImageIndex));

R_aMax = pData(I_aMax);
T_aMax = thisYearTime(I_aMax);

R_pMax = repmat(R_pMax,size(R_aMax));
T_pMax = repmat(T_pMax,size(T_aMax));

end

function pData = getCurrentPD(data0,celli,cellj,dataInfo)
pData = squeeze(data0(celli,cellj,:));
end

function cellRange = getCellRange(data0,dataInfo)
domainLen = dataInfo.arf.searchLen;
cellRange = round(size(data0,1)/2)-floor(domainLen/2):...
    round(size(data0,1)/2)+floor(domainLen/2);
end

function data0 = processData(data,duration)
% data [loc,loc,time]
% unit: [mm/h];
B = 1/duration*ones(1,1,duration);
data0 = convn(data,B,'same');
end

function dataInfo = getDataInfo(regionName,regionYear)

dataInfo = struct;
dataInfo.Years = regionYear;
dataInfo.region = getfield(REGIONS_info(),regionName);
dataInfo.ensNo = getEnsNos();
if strcmp(regionYear,'2060-2080')
    dataInfo.fileGetPath = ['D:/UKCP18_Future_2060_2080/',dataInfo.region.Name];
    dataInfo.savePath = 'D:/UKCP18/fixedAreaARF';
elseif strcmp(regionYear,'1980-2000')
    dataInfo.fileGetPath = ['D:/UKCP18/',dataInfo.region.Name];
    dataInfo.savePath = 'D:/UKCP18/fixedAreaARF';
elseif strcmp(regionYear,'2007-2018')
    dataInfo.fileGetPath = ['D:/UKCP18/',dataInfo.region.Name];
    dataInfo.savePath = 'D:/UKCP18/fixedAreaARF';
elseif strcmp(regionYear,'1990-2014')
    dataInfo.fileGetPath = ['K:/GEAR-1hr'];
    dataInfo.savePath = 'D:/UKCP18/fixedAreaARF';
end
dataInfo.figPath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP';
dataInfo.season = 'JJA';
dataInfo.scaleF = 32;
dataInfo.gridReso = 2.2;

arf = struct;
arf.duration =  [ 1, 3, 6,12,24];
arf.threshold = [20,25,30,35,40]./arf.duration;
arf.searchLen = [49,49];
arf.scanSize = [1,1];%[45,45];%[21,21];
arf.gridSize = arf.searchLen+arf.scanSize*2-2;
arf.area = [10,30,80,150,250,380,530,700,900,1250,2000,4000,8000]./(dataInfo.gridReso).^2;
dataInfo.arf = arf;

end

function [data,time] = getData(dataInfo,ensNo,yearRange)

[data,time] = deal([]);
if isempty(yearRange)
yearRange = getYearRange(dataInfo);
end
if strcmp(dataInfo.season,'JJA')
    for mon = 6:8
        if strcmp(ensNo,'RAD')
            D = load(sprintf('%s%sRadar_mon%02d.mat',dataInfo.fileGetPath,filesep,mon),...
                'Rain');
            RainEnsembles = {D.Rain};clear D;
            shapeSize = [size(RainEnsembles{1},[1,2]),eomday(1999,mon)*24,...
                size(RainEnsembles{1},3)/eomday(1999,mon)/24];
            data = cat(3,data,reshape(RainEnsembles{1},shapeSize));
            if mod(size(data,3),24*30)==720
                data(:,:,end+1:24*30,:) = []; % # need to check/
            elseif mod(size(data,3),24*30)~=0
                data(:,:,floor(size(data,3)/24/30)*24*30+1:end,:) = [];
            end
        elseif strcmp(ensNo,'GEAR')
            A = getGearMonth(dataInfo,mon,yearRange);
            shapeSize = [size(A,[1,2]),eomday(1999,mon)*24,...
                size(A,3)/eomday(1999,mon)/24];
            data = cat(3,data,reshape(A,shapeSize));
            if mod(size(data,3),24*30)==720
                data(:,:,end+1:24*30,:) = []; % # need to check/
            elseif mod(size(data,3),24*30)~=0
                data(:,:,floor(size(data,3)/24/30)*24*30+1:end,:) = [];
            end
        else
            load(sprintf('%s%sEnsems_%s_mon%02d.mat',dataInfo.fileGetPath,filesep,ensNo,mon),...
                'RainEnsembles');
            shapeSize = [size(RainEnsembles{1},[1,2]),30*24,size(RainEnsembles{1},3)/30/24];
            data = cat(3,data,reshape(RainEnsembles{1},shapeSize));
        end
    end
    data = reshape(data,[size(data,[1,2]),prod(size(data,[3,4]))]);
    
    [hh,dd,yy] = meshgrid(1:24,1:30,yearRange);
    for mon = 6:8
        T0 = datetime(yy,mon,dd,hh,0,0);
        T0 = permute(T0,[2,1,3]);
        T0 = reshape(T0(:),30*24,numel(T0)/30/24);
        time = cat(1,time,T0);
    end
    time = time(:);
    
else
    error('Check dataInfo.season');
end


% clip to dataInfo.gridsize
offset = round(dataInfo.arf.gridSize(1)/2);
rowMid = round(size(data,1)/2);
colMid = round(size(data,2)/2);
data = data(rowMid-offset:rowMid+dataInfo.arf.gridSize-offset-1,...
    colMid-offset:colMid+dataInfo.arf.gridSize-offset-1,:);

% convert to double
data = double(data)/dataInfo.scaleF;
end

function data = getGearMonth(dataInfo,mon,YEARVEC)
data = [];
region = dataInfo.region;
if isempty(YEARVEC)
    YEARVEC = 1990:2014;
end
filePath = dataInfo.fileGetPath;%'K:/GEAR-1hr';
ncFileName = [filePath,filesep,'CEH-GEAR-1hr_199501.nc'];
E = ncread(ncFileName,'x'); N = ncread(ncFileName,'y');%max to min
[N,E] = meshgrid(N/1000,E/1000);

[region.i,region.j] = getRegionIJ(E,N,region.minE,region.minN);
for year = YEARVEC
    % note: save format of gear is different from other
    % dataset, where
    % 1. northing vec is from max to min.
    % 2. 1km resolution
    [filePath,fileName] = getGear1hr_FileName(year,mon);
    listaRain = [filePath,fileName];
    rain_temp = squeeze(ncread(listaRain,'rainfall_amount',...
        [region.i,round(region.j+(-region.dimN)*region.dx+1),1],...
        [round(region.dimE*region.dx),round(region.dimN*region.dx),Inf]));
    rain_temp = rain_temp(:,end:-1:1,:);
    % # aggregate to 2.2
    % rain_temp = imresize3(rain_temp,'Scale',[5,5,1],'Method','nearest');
    % rain_temp = imresize3(rain_temp,'Scale',[1/11,1/11,1],'Method','box');
    % quicker aggregate % any how better to use the upper one if PC storage
    % allowes u to do that!
    rain_temp = imresize3(rain_temp,'Scale',[1/2.2,1/2.2,1],'Method','box');
    data = cat(3,data,rain_temp);
end
    function [filePath,fileName] = getGear1hr_FileName(year,mon)
        filePath = 'K:/GEAR-1hr/';
        fileName = ['CEH-GEAR-1hr_',sprintf('%04d%02d.nc',year,mon)];
    end

end
function [yearRange] = getYearRange(dataInfo)
% only valid for JJA.
% please check it if other season is used.
switch(dataInfo.Years)
    case '2060-2080' % cpm2.2
        yearRange = 2061:2080;
    case '2020-2040' % cpm2.2
        yearRange = 2021:2040;
    case '1980-2000' % cpm2,2
        yearRange = 1981:2000;
    case '2007-2018' % nimrod
        yearRange = 2007:2018;   
    case '1990-2014' % gear 1hr
        yearRange = 1990:2014;
end
end
