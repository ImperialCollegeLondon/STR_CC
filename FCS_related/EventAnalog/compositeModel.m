% FINAL VERSION:
clear;clc

for configVec = [0 1 0 0;
        0 0 1 0;
        0 0 0 1]
    STATS = table();
    for expNo = {'C1P5-Rad1','C1P5In-Rad1'}%{'C1P5In-Rad1','C1P5-Rad1','C1P5In-1','C1P5-1','C2P5-1','C1P5-2','C1U5-1','C1U10-1'}
        expNo = string(expNo);
        
        [version,testConfig,filefolder] = getExpNoInfo(expNo);
        testConfig.IDW = logical(configVec(1));% false is better (at least for event 90)
        testConfig.uniqueEvent = logical(configVec(2));
        testConfig.includePrevious = logical(configVec(3));
        
        savePath = 'K:\DATA_FCS\CrossVal\';
        A = load([savePath,sprintf('SummarySTATS_%s_IDW%g_uniqueEvent%g_includePrevious%g.mat',...
            expNo,testConfig.IDW,testConfig.uniqueEvent,testConfig.includePrevious)],...
            'STATS','testConfig');
        if strcmp(expNo,'C1P5In-Rad1')
            STATS = [STATS;A.STATS(A.STATS.eventTi./A.STATS.eventTT <= 0.35, :)];
        elseif strcmp(expNo,'C1P5-Rad1')
            STATS = [STATS;A.STATS(A.STATS.eventTi./A.STATS.eventTT > 0.35, :)];
        else
            error('Check Input')
        end
    end
    save([savePath,sprintf('SummarySTATS_%s_IDW%g_uniqueEvent%g_includePrevious%g.mat',...
            'C1P5RadComb',testConfig.IDW,testConfig.uniqueEvent,testConfig.includePrevious)],...
            'STATS','testConfig');
end

