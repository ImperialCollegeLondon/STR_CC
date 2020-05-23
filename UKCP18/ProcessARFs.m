%
% This script is to calculate ARFs using a storm-centered approach.
%
% Data: 1980-2000 UKCP CPM rainfall or 2060-2080 UKCP CPM rainfall
%
%
%
% ref: https://www.researchgate.net/publication/229038327_Estimation_of_
%      Average_Rainfall_Areal_Reduction_Factors_in_Texas_Using_NEXRAD_Data
%

clear;clc

% for ensNo = dataInfo.ensNo(:)
for dataInfo = [getDataInfo('SCO','1980-2000')]%% getDataInfo('SCO','2060-2080'),
    for ensNo = dataInfo.ensNo(1,3:end)%dataInfo.ensNo(:)
        ensNo = ensNo{1};
        [data,time] = getData(dataInfo,ensNo);
        [T,SNAPSHOTS,TIMES] = deal([]);
        for year = reshape(unique(time.Year),1,[])
            stats = []; % <table>
            thisYearIndex = time.Year == year;
            thisYearTime = time(thisYearIndex);
            
            for duri = 1:length(dataInfo.arf.duration)
                duration = dataInfo.arf.duration(duri);
                
                startT = tic;
                
                threshold = getThreshold(duration,dataInfo);
                data0 = processData(data(:,:,thisYearIndex),duration);% process data
                [allSnaps,allTimes] = deal([]);
                cellRange = getCellRange(data0,dataInfo);
                for celli = cellRange
                    for cellj = cellRange
                        if ~ isscalar(celli)
                            error('check for loop for celli & cellj');
                        end
                        curData = getCurrentDomain(data0,celli,cellj,dataInfo);
                        [oneSnap,oneTime] = getMaxima(curData,thisYearTime,threshold);
                        % save this snap to allSnaps
                        if ~isempty(oneSnap)
                            allSnaps{length(allSnaps)+1} = reshape(oneSnap,[size(oneSnap),1]);
                            allTimes{length(allTimes)+1} = oneTime;
                        end
                    end
                end
                %
                if ~isempty(allSnaps)
                    allSnaps = cellfun(@(x)reshape(x,[1,size(x)]),allSnaps,'UniformOutput',false);
                    allSnaps = permute(cell2mat(cat(1,allSnaps(:))),[2,3,1]);
                    allTimes = cat(1,allTimes{:});
                    
                    [allSnaps,allTimes] = filterOutSameStorm(allSnaps,allTimes);
                    
                    am = squeeze(allSnaps(round(size(allSnaps,2)/2),round(size(allSnaps,2)/2),:));
                    arf = getArf(allSnaps,allTimes,dataInfo);
                    dur = repmat(duration,size(allTimes));
                    ens = repmat(string(ensNo),size(allTimes));
                    T = [T;table(allTimes,dur,arf,ens,am)];
                    SNAPSHOTS = cat(2,SNAPSHOTS,allSnaps(:)');
                    TIMES = cat(2,TIMES,allTimes(:)');
                    fprintf('ensNo %s duration %02d %03.1f seconds \n',ensNo,duration,toc(startT));
                else
                end
            end
            
        end
        SNAPSHOTS = SNAPSHOTS';TIMES = TIMES';
        save(sprintf('%s%sARFs_%s_ensNo%s_%s.mat',dataInfo.savePath,filesep,...
            dataInfo.region.Name,ensNo,dataInfo.Years),'T','TIMES','SNAPSHOTS','-v7.3');
    end
    
end
%% AUXILLARY FUNCTION
function thisarf = getArf(allSnaps,allTimes,dataInfo)
thisarf = [];
for ii = 1:size(allSnaps,3)
    thisImage = squeeze(allSnaps(:,:,ii));
    tag = 0;
    % regionprops('table',logical(thisImage),'MajorAxisLength','MinorAxisLength')
    for area = reshape(dataInfo.arf.area,1,[])
        tag = tag+1;
        [centralI,centralJ,thisHighImage,area0] = processImage(thisImage,area);
        fitRes = fitIt(thisHighImage,area0,centralI,centralJ);
        thisarf(ii,tag) = fitRes(3)./nanmax(thisHighImage(:));  
    end
    thisarf(ii,tag);
    
    hold on;
    thisImage0 = thisImage;thisImage0(thisImage0 == 0) = NaN;
    pcolor(1:size(thisImage,1),1:size(thisImage,2),thisImage0);shading flat;
    plot(size(thisImage,1)/2+1,size(thisImage,1)/2+1,'rx');
    bb = sqrt(area/pi/fitRes(2));aa = fitRes(2)*bb;
    plot([size(thisImage,1)/2+1,size(thisImage,1)/2+1+aa*cos(fitRes(1)),...
        size(thisImage,1)/2+1+bb*sin(fitRes(1)),size(thisImage,1)/2+1],...
        [size(thisImage,1)/2+1,size(thisImage,1)/2+1+aa*sin(fitRes(1)),...
        size(thisImage,1)/2+1-bb*cos(fitRes(1)),size(thisImage,1)/2+1],...
        'r-','linewidth',2);
    xlim([1,size(thisImage,1)]);ylim([1,size(thisImage,2)]);
    cptcmap('cw1-013','mapping','scaled','ncol',10,'flip',true);
    caxis([0,30]);
    hold off;
    close all
    pause(0.01)
        
%     pause(0.1)
end
    function [centralI,centralJ,thisHighImage,area0] = processImage(thisImage,area)
        if area < 800
            scaleF = 2.2;% 2.2;%2.2;
            thisHighImage = imresize(thisImage,scaleF,'box');
            area0 = area*scaleF^2;
        else
            thisHighImage = thisImage;
            area0 = area;
        end
        centralI = round(size(thisHighImage,1)/2);
        centralJ = round(size(thisHighImage,2)/2);
    end
    function fitRes = fitIt(thisHighImage,area,centralI,centralJ)
        % [theta,aspect,fval]
        % xi --> rowi;  yi: --> coli;
        % theta: angle to y-axis (clockwise)
        % included points need to meet:
        % yi = (coli-centralJ);
        % xi = (rowi-centralI);
        % (yi*cos(theta)+xi*sin(theta)).^2./(a^2)+(yi*cos(theta)-xi*sin(theta)).^2./(b^2)<1
        % at the mean time: area = pi*a*b --> b = sqrt(area/pi/R);a = R*b;
        % Boudary: theta = [0,pi]; R = [1,10];
        % optimized function:
        % to max(obs(inElip));
        
%         if area < 20
%             lb = [0,1];
%             ub = [pi,2]; % 'ConstraintTolerance',0.05
%         else
%             lb = [0,1];
%             ub = [pi,10]; % 'ConstraintTolerance',0.05
%         end
lb = [0,1];
ub = [pi,9]; % 'ConstraintTolerance',0.05

        options = optimoptions('ga','Display','off','TimeLimit',240*60,...
            'HybridFcn',@fmincon,'PopulationSize',200,'MaxGenerations',200);
        % gaoptimset();
        [yi,xi] = meshgrid(((1:size(thisHighImage,2)))-(size(thisHighImage,1)/2),...
            (1:size(thisHighImage,1))-(size(thisHighImage,1)/2));
        [res,fval,exitflag,output,population,scores] = ga(@(X)objFunc(X(1),X(2),area,yi,xi), ...
            2, [],[],[],[],lb,ub,[],options);
        
        theta = res(1);
        R = res(2);
        b = sqrt(area/pi/R);a = R*b;
        inElip = (yi(:)*cos(theta)+xi(:)*sin(theta)).^2./(a^2)+...
            (yi(:)*sin(theta)-xi(:)*cos(theta)).^2./(b^2)<1;
        if area>0%25
            thisAreaRain = -fval/area;
        else
            thisAreaRain = -fval/nansum(inElip(:));
        end
        fitRes = [res,thisAreaRain];
    end
    function res = objFunc(theta,R,area,yi,xi)
        b = sqrt(area/pi/R);a = R*b;
        inElip = (yi(:)*cos(theta)+xi(:)*sin(theta)).^2./(a^2)+...
            (yi(:)*sin(theta)-xi(:)*cos(theta)).^2./(b^2)<1;
        res = nansum(-1*thisHighImage(inElip(:)));%./nansum(inElip(:));
    end
end
function threshold = getThreshold(duration,dataInfo)
threshold = dataInfo.arf.threshold(dataInfo.arf.duration == duration);
end

function [allSnaps,allTimes] = filterOutSameStorm(allSnaps,allTimes)
centralTS = getCentral(allSnaps);
[~,I] = sort(centralTS,'descend');
allTimes = allTimes(I);
allSnaps = allSnaps(:,:,I);
[~,IA,~] = unique(round(datenum(allTimes)*(24/24)),'stable');
allTimes = allTimes(IA);
allSnaps = allSnaps(:,:,IA);
    function centralTS = getCentral(curData)
        ci = round(size(curData,2)/2);
        centralTS = squeeze(curData(ci,ci,:));
    end
end

function [oneSnap,oneTime] = getMaxima(curData,thisYearTime,threshold)
% get maximum at the premise that the storm centred here
centralTS = getCentral(curData);
[~,I] = sort(centralTS,'descend');
okFlag = false;
thisImageIndex = 0;
if okFlag == false && thisImageIndex < size(curData,3)
    thisImageIndex = thisImageIndex+1;
    oneSnap = squeeze(curData(:,:,I(thisImageIndex)));
    if centralTS(I(thisImageIndex)) >= threshold && nanmax(oneSnap(:)) == centralTS(I(thisImageIndex))
        [i_temp,j_temp] = find(oneSnap == nanmax(oneSnap(:)),1);
        if i_temp == j_temp && i_temp == round(size(curData,2)/2)
            okFlag = true;
        end
    end
end
if okFlag
    thisImageIndex;
    % stats = regionprops('table',oneSnap,'Area','MaxIntensity');
    oneTime = thisYearTime(I(thisImageIndex));
else
    [oneSnap,oneTime] = deal([]);
end
    function centralTS = getCentral(curData)
        ci = round(size(curData,2)/2);
        centralTS = squeeze(curData(ci,ci,:));
    end
end

function curData = getCurrentDomain(data0,celli,cellj,dataInfo)
offset = floor(dataInfo.arf.scanSize/2);
rowIndex = celli-offset:celli+(dataInfo.arf.scanSize-offset)-1;
colIndex = cellj-offset:cellj+(dataInfo.arf.scanSize-offset)-1;
curData = data0(rowIndex,colIndex,:);
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
    dataInfo.savePath = 'D:/UKCP18/stormcentredARF';
elseif strcmp(regionYear,'1980-2000')
    dataInfo.fileGetPath = ['D:/UKCP18/',dataInfo.region.Name];
    dataInfo.savePath = 'D:/UKCP18/stormcentredARF';
end
dataInfo.figPath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP';
dataInfo.season = 'JJA';
dataInfo.scaleF = 32;
dataInfo.gridReso = 2.2;

arf = struct;
arf.duration =  [ 1, 3, 6,12,24];
arf.threshold = [20,25,30,35,40]./arf.duration;
arf.searchLen = [49,49];%;[101,101];%[51,51];
arf.scanSize = [45,45];%[21,21];
arf.gridSize = arf.searchLen+arf.scanSize*2-2;
arf.area = [10,30,80,150,250,380,530,700,900]./(dataInfo.gridReso).^2;
dataInfo.arf = arf;

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

function [yearRange] = getYearRange(dataInfo)
% only valid for JJA.
% please check it if other season is used.
yearRange = strcmp(dataInfo.Years,'2060-2080')*(2061:2080)+...
    strcmp(dataInfo.Years,'2020-2040')*(2021:2040)+...
    strcmp(dataInfo.Years,'1980-2000')*(1981:2000);
end




