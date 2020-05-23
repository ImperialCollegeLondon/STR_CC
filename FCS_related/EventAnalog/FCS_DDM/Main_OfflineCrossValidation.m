function STATS = Main_OfflineCrossValidation(version,testConfig,filefolder)
%
% Main_OfflineCrossValidation will evluate the statistics related to
% testConfig using the data saved in filefolder.
% 
% Example:
% [version,testConfig,filefolder] = getExpNoInfo('C1U5-1');
% STATS = Main_OfflineCrossValidation(version,testConfig,filefolder);
%
%
% @ Yuting Chen
% yuting.chen17@imperial.ac.uk


mkdir(filefolder);

for eventNo4val = 1:157
    % Save result
        try
            [LIB] = Main_OfflineTraining(version,eventNo4val,testConfig);
            [STATS,FOREOUTPUT] = Main_OfflineTesting(LIB,version,testConfig);
            save([filefolder,'STATS_',sprintf('%s_ev%03d.mat',version,eventNo4val)],...
                'STATS','FOREOUTPUT','version');
            aux_modelSetting(FOREOUTPUT,STATS)
        catch me
            me
            eventNo4val
        end
end

end



