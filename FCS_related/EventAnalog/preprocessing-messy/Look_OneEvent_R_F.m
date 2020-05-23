% --------------------------------------------------------- %
% See One event
% Both Rainfall and Flood Map
%
% Note:
% It takes a quite while to open each Ascii file (~200Mb)
% 52.4 GB * 4 in total
%
% @ Yuting Chen
% Update: 2020.01.09
% --------------------------------------------------------- %

% clear;clc
% close all

%% Configuration
FloodPath = 'G:\BIGDATA\TOPIC 2\FileFromLIPEN\Flood Maps\Batch_1\';
RainfallPath = 'G:\BIGDATA\TOPIC 2\FileFromLIPEN\Rainfall_in_157Floods\KED_ASCII\';
EventNo = '037';
EventTimes = '201408011300-201408021535';


%% Get Info of all Files
FloodFNs = getFloodFiles(FloodPath,EventNo);
RainfallFNs = getRainfallFiles(RainfallPath,EventTimes);

% Import Files
try
    [FMap,FTag] = ImportRPSFloodMap(FloodFNs);
    [RMap,RTag] = ImportRainfallMap(RainfallFNs);
catch
    1;
end

save(['Event',EventNo,'.mat'], 'FMap', 'RMap', 'FTag', 'RTag');

%% Show Rainfall
load(['Event',EventNo,'.mat'], 'FMap', 'RMap', 'FTag', 'RTag');
R = zeros(34,30);

for wc = 1:length(RMap)
    
    % R = R+RMap{wc};
    pcolor(RMap{wc});
    shading flat;
    cptcmap('GMT_drywet', 'mapping', 'direct');
    colorbar;
    pause(0.02);
    
end

%% Show FloodMap
for wc = 1:length(FMap)
    
    pcolor(exp(FMap{wc})-1);
    shading flat;
    cptcmap('GMT_drywet', 'mapping', 'scaled');
    colorbar;
    pause(0.5);%
    
end

%% AUXILLARY FUNCTION

function FileNames = getFloodFiles(Path,TagNo)

D = dir([Path,TagNo,'_WC*.asc']);
FileNames = arrayfun(@(x)strcat(string(x.folder),'\',string(x.name)), ...
    D, 'UniformOutput', false);

end

function FileNames = getRainfallFiles(Path,EventTimes)

startTime = datetime(EventTimes(1:12),'InputFormat','yyyyMMddHHmm');
endTime = datetime(EventTimes(14:end),'InputFormat','yyyyMMddHHmm');
FileNames = [];
tag = 1;
for No = startTime:1/24/12:endTime
    FileNames{tag} = strcat(Path,string(No,'yyyyMMddHHmm'),'.ASC');
    tag = tag+1;
end

end


function [FMap,FTag] = ImportRPSFloodMap(FloodFNs)

try
    FMap = [];
    for Num = 1:numel(FloodFNs)
        [OUT, ncols, nrows, xllcorner, yllcorner, cellsize, nodata] = ...
            ascii_reader(FloodFNs{Num});
        FMap{Num} = imresize(OUT,1/25);
        FTag{Num} = FloodFNs{Num};
    end
catch
    1;
end

end

function [RMap,RTag] = ImportRainfallMap(RainfallFNs)

try
    RMap = [];
    for Num = 1:numel(RainfallFNs)
        [OUT, ncols, nrows, xllcorner, yllcorner, cellsize, nodata] = ...
            ascii_reader(RainfallFNs{Num});
        RMap{Num} = OUT;
        RTag{Num} = RainfallFNs{Num};
    end
catch
    1;
end

end




