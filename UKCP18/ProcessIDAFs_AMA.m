%
% This script is to derive required data for IDAF based on annual maxima ananlysis
%
% Data: 1980-2000 UKCP CPM rainfall or 2060-2080 UKCP CPM rainfall
%
%
% ref: ...2002
% based on C. De Michele 2011 func.6 in Journal of Hydrology
%
% @ Yuting CHEN
% yuting.chen17@imperial.ac.uk
% update: 2020/08/06

clear;clc
close all

dataInfo = getDataInfo('CPM_NE','1980-2000');
for ensNo = dataInfo.ensNo(1:end)
    ensNo = ensNo{1};
    for dataInfo = [
            getDataInfo('CPM_NE','1980-2000'),getDataInfo('CPM_NE','2060-2080'),...
            getDataInfo('CPM_NW','1980-2000'),getDataInfo('CPM_NW','2060-2080')]
        % getDataInfo('CPM_S','1980-2000'),getDataInfo('CPM_S','2060-2080')]
        % [getDataInfo('EUK','1980-2000'),getDataInfo('EUK','2060-2080')]
        % [getDataInfo('SCO','1980-2000'),getDataInfo('SCO','2060-2080'),...
        %    getDataInfo('WAL','1980-2000'),getDataInfo('WAL','2060-2080'),...
        %     getDataInfo('EUK','1980-2000'),getDataInfo('EUK','2060-2080')]
        
        saveFile = sprintf('%s%sIDAFs_%s_ensNo%s_%s.mat',dataInfo.savePath,filesep,...
            dataInfo.region.Name,ensNo,dataInfo.Years);
        if fopen(saveFile)~=-1
            % fclose('all')
            continue;
        end
        
        [T,ALLRMAX,TIMES] = deal([]);
        
        [~,time] = getData(dataInfo,ensNo,'noData');% #to do# func can be modifed for quicker reading
        [data,time] = getData(dataInfo,ensNo,'');
        for year = reshape(unique(time.Year),1,[])
            
            stats = []; % <table>
            
            tic
            thisYearIndex = time.Year == year;
            thisYearTime = time(thisYearIndex);
            dataThisYear = data(:,:,thisYearIndex);
            % clear data;
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
                    if ~isempty(allT_Max)
                        T = [T;table(allT_Max,dur,ens,am,area)];
                        TIMES = cat(2,TIMES,allT_Max(:)');
                    end
                    fprintf('ensNo%s year%04d duration %02d area %04.1d %02.1f seconds \n',...
                        ensNo,year,duration,area,toc(startT));
                end
            end
            
        end
        ALLRMAX = ALLRMAX';TIMES = TIMES';
        save(saveFile,'T','TIMES','-v7.3');
        clear data dataThisYear
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
idaf.area = [1:42].^2*(dataInfo.gridReso).^2;
dataInfo.idaf = idaf;

end

function [data,time] = getData(dataInfo,ensNo,noData)
% has not covert to double
[data,time] = deal([]);
yearRange = getYearRange(dataInfo);
if strcmp(dataInfo.season,'JJA')
    
    [hh,dd,yy] = meshgrid(1:24,1:30,yearRange);
    for mon = 6:8
        T0 = datetime(yy,mon,dd,hh,0,0);
        T0 = permute(T0,[2,1,3]);
        T0 = reshape(T0(:),30*24,numel(T0)/30/24);
        time = cat(1,time,T0);
    end
    time = time(:);
    
    if ~strcmp('noData',noData)
        for mon = 6:8
            load(sprintf('%s%sEnsems_%s_mon%02d.mat',dataInfo.fileGetPath,filesep,ensNo,mon),...
                'RainEnsembles');
            shapeSize = [size(RainEnsembles{1},[1,2]),30*24,size(RainEnsembles{1},3)/30/24];
            data = cat(3,data,reshape(RainEnsembles{1},shapeSize));
        end
        data = reshape(data,[size(data,[1,2]),prod(size(data,[3,4]))]);
        clear RainEnsembles
    end
    % trim those points outside UK
    tic
    UKMap = load('H:\CODE_MATLAB\SpatialTemporalDATA\shapeFileFolder\UKBorderGrid.mat');
    [E1,N1] = getEN(dataInfo.region);
    [in] = inpolygon(E1,N1,UKMap.borderE/1000,UKMap.borderN/1000);
    
    data = single(data)/dataInfo.scaleF;% convert to double here.
    
    in = repmat(in,[1,1,size(data,3)]);
    data(~in) = NaN;
    toc
    
else
    error('Check dataInfo.season');
end
% if ~strcmp('noData',noData) & ~isRegularRegion(dataInfo.region)
%     data = padarray(data,[1 1,0],NaN,'both');
%     N1 = padarray(N1,[1 1],NaN,'both');
%     E1 = padarray(E1,[1 1],NaN,'both');
%     angle = atan((dataInfo.region.E(2)-dataInfo.region.E(1))/...
%         (dataInfo.region.N(2)-dataInfo.region.N(1)))*180/pi;
%     E1 = imrotate(E1,-angle,'nearest','loose');
%     N1 = imrotate(N1,-angle,'nearest','loose');
%     data = imrotate(data,-angle,'nearest','loose');
% end
end

function [yearRange] = getYearRange(dataInfo)
% only valid for JJA.
% please check it if other season is used.
yearRange = strcmp(dataInfo.Years,'2060-2080')*(2061:2080)+...
    strcmp(dataInfo.Years,'2020-2040')*(2021:2040)+...
    strcmp(dataInfo.Years,'1980-2000')*(1981:2000);
end
