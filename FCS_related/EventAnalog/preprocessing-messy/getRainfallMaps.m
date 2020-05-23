% ----------------------------------------------------------------------- %
% Ger Rain maps for each rainfall event
% 
% 
% Yuting Chen
% Update: 2020.01.31
% ----------------------------------------------------------------------- %


%% Configuration

REGIONS = REGIONS_info();
region = REGIONS.Birm;

% FloodPath = 'G:\BIGDATA\TOPIC 2\FileFromLIPEN\Flood Maps\Batch_4\';
RainfallPath = 'G:\BIGDATA\TOPIC 2\FileFromLIPEN\Rainfall_in_157Floods\KED_ASCII\';


%% Get and Merge Info of FloodMap Files for this EventNo

D = dir([FloodPath,'*_WC*_PEAK.asc']);

EventNames = unique(cellstr(arrayfun(@(x)x.name(1:3), ...
    D, 'UniformOutput', false)));

for EventNo = EventNames'
    
    EventNo = EventNo{1};
    load(sprintf('G:/BIGDATA/TOPIC 2/ProcessedFiles/FloodMaps/FloodMaps_Merged_EventNo%s.mat',...
        EventNo),'FMap','FMapFName','FMapInfo','E','N','unit');
    
end

%% Get and Merge Info of FloodMap Files for this EventNo

for no = 1:157
    try
        
        EventNo = sprintf('%03d',no);
        EventTimes = getEventTimes(EventNo);%'201408011300-201408021535';
        RainfallFNs = getRainfallFiles(RainfallPath,EventTimes);
        [RMap,E,N,RTag] = ImportRainfallMap(RainfallFNs);
        [EE,NN] = meshgrid(E/1000,N/1000);
        % RMap = cellfun(@(x)imresize(x,1000/10,'nearest'),RMap,'UniformOutput',false);
        
        save(sprintf('G:/BIGDATA/TOPIC 2/ProcessedFiles/RainEvent/EventNo%s.mat',EventNo),...
            'RMap','EE','NN','RTag','EventTimes');
        fprintf('%s: Done \n',EventNo);
        
        
    catch me
        fprintf('No Event Info this time \n');
    end
    
end




%% Show Rainfall
% load(['Event',EventNo,'.mat'], 'FMap', 'RMap', 'FTag', 'RTag');
R = zeros(34,30);

for wc = 1:length(RMap)
    
    % R = R+RMap{wc};
    RMap{wc}((RMap{wc}==0)) = NaN;
    pcolor(RMap{wc});
    shading flat;
    cptcmap('precip_meteoswiss', 'mapping', 'scaled');
    caxis([0,20])
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

function EventTimes = getEventTimes(EventNo)
% EventNo = '037';
% EventTimes = getEventTimes(EventNo);%'201408011300-201408021535';

eventNo = str2num(EventNo);
filePath = 'G:\BIGDATA\TOPIC 2\FileFromLIPEN\Rainfall_in_157Floods\';
[~, text, ~] = xlsread([filePath,'Birmingham_SelectedEvents_DataDriven.csv'], 1, ...
    sprintf('F%d:F%d',eventNo+1,eventNo+1));
EventTimes = text{1}(end-28:end-4);

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


function [RMap,E,N,RTag] = ImportRainfallMap(RainfallFNs)

try
    RMap = zeros(34,30,0);
    for Num = 1:numel(RainfallFNs)
        [OUT, ncols, nrows, xllcorner, yllcorner, cellsize, nodata] = ...
            ascii_reader(RainfallFNs{Num});
        RMap = cat(3,RMap,OUT);
        RTag{Num} = RainfallFNs{Num};
    end
    E = xllcorner:cellsize:xllcorner+cellsize*(ncols-1);
    N = yllcorner:cellsize:yllcorner+cellsize*(nrows-1);
catch
    1;
end

end








