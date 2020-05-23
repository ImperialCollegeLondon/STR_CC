function [startTime,endTime,eventNos] = getHistEventTime(source)
% getHistEventTime(source) ....
%
% Example: 
%          getHistEventTime('RAD') 
% @Yt

arguments
    source {mustBeMember(source,{'KED','RAD'})}
end

switch(source)
    case 'KED'
        filePath = 'K:\DATA_FCS\RainEvent\';
        [~, text, ~] = xlsread([filePath,'Birmingham_SelectedEvents_DataDriven.csv']);
        
        EventTimes = cellfun(@(x)x(end-28:end-4),text(2:end,6),'UniformOutput', false);
        startTime = cell2mat(cellfun(@(x)datenum(datetime(x(1:12),'InputFormat','yyyyMMddHHmm')),EventTimes,...
            'UniformOutput', false));
        endTime = cell2mat(cellfun(@(x)datenum(datetime(x(14:end),'InputFormat','yyyyMMddHHmm')),EventTimes,...
            'UniformOutput', false));
        eventNos = [1:length(endTime)]';
        
    case 'RAD'
        A = readtable('K:\DATA_FCS\RainEvent\Birmingham_SelectedEvents_DataDriven_RadarVersion.csv');
        startTime = string(A.startTime);
        startTime = cell2mat(cellfun(@(x)datenum(datetime(x,'InputFormat','yyyyMMddHHmm')),startTime,...
            'UniformOutput', false));
        endTime = string(A.endTime);
        endTime = cell2mat(cellfun(@(x)datenum(datetime(x,'InputFormat','yyyyMMddHHmm')),endTime,...
            'UniformOutput', false));
        eventNos = A.eventNos;
        
    otherwise
        
end
end