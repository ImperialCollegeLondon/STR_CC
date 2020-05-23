function [STATS,FOREOUTPUT] = Main_OfflineTesting(LIB_ori,version,testConfig)

tic

[LIB] = splitDatasets(LIB_ori);


[AnalogInput] = getAnalogInput(LIB,version,testConfig);


FANum = 120;
[RES_analog,DIS_analogL1,DIS_analogL2] = Model_Analogue(AnalogInput,FANum);


% FORECASTING RAINFALL EVENT
[FOR_test,FOR_analog] = Model_ensembleForecasting(RES_analog,AnalogInput,LIB);


% DERIVING ALL EVENT-BASED FLOODMAPS AND OBSERVED FLOOD MAPS
[FLO_test,FLO_pred,E,N] = Model_floodPrediction(FOR_test,FOR_analog);

toc

%-------------- CONSTRUCTING OUTPUT --------------%

% [STATS] = evaluateFloodPrediction(FLO_test,FLO_pred,E,N,DIS_analogL1,testConfig);
[STATS] = evaluateFloodPrediction_FCSversion(FLO_test,FLO_pred,E,N,DIS_analogL1,testConfig);

% ORIGINAL FORECAST OUTPUT
[FOREOUTPUT] = formatOutput(AnalogInput,RES_analog,DIS_analogL1,DIS_analogL2,FOR_test,FOR_analog,FLO_test,FLO_pred,E,N);

end


function [LIB] = splitDatasets(LIB)
LIB = LIB;
end

function [AnalogInput] = getAnalogInput(LIB,version,testConfig)
%% configure input

AnalogInput = struct;

AnalogInput.RIMAGES = LIB.RainLib.train.RIMAGES;
AnalogInput.CRIMAGES = LIB.RainLib.train.CRIMAGES;
AnalogInput.RIMAGESval = LIB.RainLib.test.RIMAGES;
AnalogInput.CRIMAGESval = LIB.RainLib.test.CRIMAGES;
AnalogInput.raintime = LIB.RainLib.train.raintime;
AnalogInput.raintimeval = LIB.RainLib.test.raintime;

try AnalogInput.timegap = testConfig.analogtimegap;catch
AnalogInput.timegap = 0;
end

updatePredClimate(version,testConfig);
updatePredRain(version,testConfig);


    function updatePredRain(version,testConfig)
        RIMAGES = AnalogInput.RIMAGES;
        RIMAGESval = AnalogInput.RIMAGESval;
        
        PREDRAIN = LIB.RainLib.train.PREDRAIN;
        PREDRAINval = LIB.RainLib.test.PREDRAIN;
        
        
        [PREDRAIN_pca_lm,PREDRAINval_pca_lm,U_pm,~] = getImageFeature(RIMAGES,RIMAGESval,...
            testConfig.imageMethod,testConfig.pcNum);
        
  
        % XTrain = transpose(preTrans(CRIMAGES));
        % % Preprocessing of PCA func:
        % % original ingredients data centered by subtracting the column means from corresponding columns.
        % [coeff,scoreTrain,latent,tsquared,explained,mu]  = pca(XTrain,'NumComponents',5);
        % U_pm = coeff(:,1:1);% 95% variance is preserved
        % X_lp = exp(transpose(preTrans(CRIMAGES)))-0.01;
        % PREDCRAIN_pca_lm = X_lp*U_pm;
        % X_lp = exp(transpose(preTrans(CRIMAGESval)))-0.01;
        % PREDCRAINval_pca_lm = X_lp*U_pm;
        if ~isfield(testConfig,'excludeWAR') || testConfig.excludeWAR == false
            AnalogInput.PREDRAIN  = [PREDRAIN_pca_lm,PREDRAIN.WAR'];%,PREDRAIN.IMF'];
            AnalogInput.PREDRAINval = [PREDRAINval_pca_lm,PREDRAINval.WAR'];%,PREDRAINval.IMF'];
        else
            AnalogInput.PREDRAIN  = [PREDRAIN_pca_lm];%PREDRAIN.WAR',PREDRAIN.IMF'];
            AnalogInput.PREDRAINval = [PREDRAINval_pca_lm];%,PREDRAINval.WAR',PREDRAINval.IMF'];
        end
    end

    function updatePredClimate(version,testConfig)
        
        switch (version)
            case 'v1'
                PRED = LIB.ClimateLib.train;
                AnalogInput.PRED = [PRED.U',PRED.V',PRED.RH',PRED.DT',PRED.GEOZ']; % PRED.MSP',%[PRED.T',PRED.CRWC'];
                
                PREDval = LIB.ClimateLib.test;
                AnalogInput.PREDval = [PREDval.U',PREDval.V',PREDval.RH',PREDval.DT',PREDval.GEOZ']; % PREDval.MSP %[PREDval.T',PREDval.CRWC'];
            case 'v2'
                PRED = LIB.ClimateLib.train;
                AnalogInput.PRED = [PRED.U',PRED.V',PRED.RH',PRED.DT',PRED.GEOZ']; % PRED.MSP',%[PRED.T',PRED.CRWC'];
                
                PREDval = LIB.ClimateLib.test;
                AnalogInput.PREDval = [PREDval.U',PREDval.V',PREDval.RH',PREDval.DT',PREDval.GEOZ']; % PREDval.MSP %[PREDval.T',PREDval.CRWC'];
            
            otherwise
        
        end
    end

end

function [FOREOUTPUT] = formatOutput(AnalogInput,RES_analog,DIS_analogL1,DIS_analogL2,FOR_test,FOR_analog,FLO_test,FLO_pred,E,N)


FOREOUTPUT = struct;
FOREOUTPUT.analogDistL1 = DIS_analogL1;
FOREOUTPUT.analogDistL2 = DIS_analogL2;
FOREOUTPUT.obsRIMAGES = AnalogInput.RIMAGESval;
FOREOUTPUT.analogRIMAGES = AnalogInput.RIMAGES(RES_analog);
FOREOUTPUT.obsTime = AnalogInput.raintimeval;
FOREOUTPUT.analogTime = AnalogInput.raintime(RES_analog);
FOREOUTPUT.obsEvNo = FLO_test.eventNo';
FOREOUTPUT.analogPredEvNo = FLO_pred.eventNo;
FOREOUTPUT.FLO_test = FLO_test;
FOREOUTPUT.FLO_pred = FLO_pred;
FOREOUTPUT.E = E;
FOREOUTPUT.N = N;

end