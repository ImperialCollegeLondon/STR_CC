function [RIMAGES,CRIMAGES,time,RIMAGESval,CRIMAGESval,timeval,EE,NN] = getRainEvents(year,eventNo,eventNosval,filePath)

% filePath = 'K:\DATA_FCS\RainEvent\Birm-34-30'; or Birm-34-30-Radar
FNS = dir([filePath,filesep,'EventNo*.mat']);
[RIMAGES,CRIMAGES,RIMAGESval,CRIMAGESval] = deal(zeros(34,30,0));
[time,timeval] = deal([]);
excluNo = eventNosval;

for fi = 1:length(FNS)
    if ismember(fi,eventNo)
        % for training
        A = load([FNS(fi).folder,filesep,FNS(fi).name],'RMap','EventTimes','EE','NN');
        RIMAGES = cat(3,RIMAGES,A.RMap);
        CRIMAGES = cat(3,CRIMAGES,cumsum(A.RMap,3));
        time = [time;getTimeVec(A.EventTimes)];
        EE = A.EE;
        NN = A.NN;
    elseif ismember(fi,excluNo)
        % for calibrating
        A = load([FNS(fi).folder,filesep,FNS(fi).name],'RMap','EventTimes');
        RIMAGESval = cat(3,RIMAGESval,A.RMap);
        CRIMAGESval = cat(3,CRIMAGESval,cumsum(A.RMap,3));
        timeval = [timeval;getTimeVec(A.EventTimes)];
    else
        errordlg('Check function:getRainEvents', 'Error Dialog');
    end
end


    function timeVec = getTimeVec(timeStr)
        % '201407081015-201407082000'
        timeVec = datetime(timeStr(1:12),'format','yyyyMMddHHmm'):minutes(5):...
            datetime(timeStr(14:end),'format','yyyyMMddHHmm');
        timeVec = reshape(timeVec,[],1);
        %
    end

end
