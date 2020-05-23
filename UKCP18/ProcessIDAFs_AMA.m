%
% This script is to derive required data for IDAF based on annual maxima ananlysis
%
% Data: 1980-2000 UKCP CPM rainfall or 2060-2080 UKCP CPM rainfall
%
%
%
% ref: ...2002
% based on C. De Michele 2011 func.6 in Journal of Hydrology
%
% @ Yuting CHEN
% yuting.chen17@imperial.ac.uk


clear;clc
close all

dataInfo = getDataInfo('WAL','1980-2000');
for ensNo = dataInfo.ensNo(1:end)
    ensNo = ensNo{1};
    for dataInfo = [getDataInfo('SCO','1980-2000'),getDataInfo('SCO','2060-2080'),...
            getDataInfo('WAL','1980-2000'),getDataInfo('WAL','2060-2080'),...
            getDataInfo('EUK','1980-2000'),getDataInfo('EUK','2060-2080')]
        
        [T,ALLRMAX,TIMES] = deal([]);
        
        [~,time] = getData(dataInfo,ensNo);% #to do# func can be modifed for quicker reading
        
        for year = reshape(unique(time.Year),1,[])

            stats = []; % <table>
            
            tic
            [data,time] = getData(dataInfo,ensNo);
            thisYearIndex = time.Year == year;
            thisYearTime = time(thisYearIndex);
            dataThisYear = data(:,:,thisYearIndex);
            clear data;
            dataThisYear = double(dataThisYear);
            toc
            
            for duri = 1:length(dataInfo.idaf.duration)
                duration = dataInfo.idaf.duration(duri);% now only look at 1 hour
                % Note: current script can be easily changed into a coarser duration
                % No need any change on data processing!
                
                for area = dataInfo.idaf.area % [1:10:100].^2*2.2*2.2;% unit:km^2
                    gridNum = round(area/2.2/2.2);% unit: -
                    
                    startT = tic;
                    
                    data0 = processData(dataThisYear,duration,gridNum);% process data
                    I = find(data0 == nanmax(data0(:)),1);
                    [loci,locj,ti] = ind2sub(size(data0),I);
                    allR_Max = nanmax(data0(:));
                    allT_Max = thisYearTime(ti);
                    allLoc_Max = NaN;%note here is not directly equal to [loci,locj]! due to valid convn;
                    
                    am = allR_Max;
                    dur = duration;
                    ens = string(ensNo);
                    T = [T;table(allT_Max,dur,ens,am,area)];
                    TIMES = cat(2,TIMES,allT_Max(:)');
                    fprintf('ensNo%s year%04d duration %02d area %04.1d %02.1f seconds \n',...
                        ensNo,year,duration,area,toc(startT));
                end
            end
            
        end
        ALLRMAX = ALLRMAX';TIMES = TIMES';
        save(sprintf('%s%sIDAFs_%s_ensNo%s_%s.mat',dataInfo.savePath,filesep,...
            dataInfo.region.Name,ensNo,dataInfo.Years),'T','TIMES','-v7.3');
    end
    
end

%% AUXILLARY FUNCTION

function data0 = processData(data,duration,gridNum)
% data [loc,loc,time]
% unit: [mm/h];
locNo = sqrt(gridNum);
if abs(locNo-round(locNo))>0.1
    error('gridNum for convolution is not interger. Please change the area size you chose!');
end
locNo = round(locNo);
B = 1/duration*ones(locNo,locNo,duration)./locNo/locNo;
% data0 = convn(data,B,'same');
data0 = convn(data,B,'valid');
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
idaf.area = [1:50].^2*(dataInfo.gridReso).^2;
dataInfo.idaf = idaf;

end

function [data,time] = getData(dataInfo,ensNo)
[data,time] = deal([]);
yearRange = getYearRange(dataInfo);
if strcmp(dataInfo.season,'JJA')
    
    for mon = 6:8
        load(sprintf('%s%sEnsems_%s_mon%02d.mat',dataInfo.fileGetPath,filesep,ensNo,mon),...
            'RainEnsembles');
        shapeSize = [size(RainEnsembles{1},[1,2]),30*24,size(RainEnsembles{1},3)/30/24];
        data = cat(3,data,reshape(RainEnsembles{1},shapeSize));
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
    
    
    % # to do #
    % # trim those points outside UK #
    tic
    UKMap = load('H:\CODE_MATLAB\SpatialTemporalDATA\shapeFileFolder\UKBorderGrid.mat');
    [E1,N1] = getEN(dataInfo.region);
    [in] = inpolygon(E1,N1,UKMap.borderE/1000,UKMap.borderN/1000);
    
    in = repmat(in,[1,1,size(data,3)]);
    data(~in) = NaN;
    toc
    
else
    error('Check dataInfo.season');
end

% convert to double
data = single(data)/dataInfo.scaleF;

end

function [yearRange] = getYearRange(dataInfo)
% only valid for JJA.
% please check it if other season is used.
yearRange = strcmp(dataInfo.Years,'2060-2080')*(2061:2080)+...
    strcmp(dataInfo.Years,'2020-2040')*(2021:2040)+...
    strcmp(dataInfo.Years,'1980-2000')*(1981:2000);
end
