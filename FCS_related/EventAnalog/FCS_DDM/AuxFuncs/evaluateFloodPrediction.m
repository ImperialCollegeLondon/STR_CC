function [STATS] = evaluateFloodPrediction(FLO_test,FLO_pred,E,N,DIS_analog,testConfig)
%
% reference: (Wing et al., 2017)
%

STATS = struct;
unit = '3km';
aggUnit = 3;
thr = 0;


D1H1 = @(sim,obs)1e-6+sum(sim(:)>thr & obs(:)>thr);
D0H1 = @(sim,obs)1e-6+sum(sim(:)<=thr & obs(:)>thr);
D1H0 = @(sim,obs)1e-6+sum(sim(:)>thr & obs(:)<=thr);
D0H0 = @(sim,obs)1e-6+sum(sim(:)<=thr & obs(:)<=thr);


RMSE = @(sim,obs)sqrt(nanmean((squeeze(sim(:))-squeeze(obs(:))).^2));
MAPE = @(sim,obs)nanmean(abs(sim(obs(:)>thr | sim(:)>thr) - ...
    obs(obs(:)>thr | sim(:)>thr)./obs(obs(:)>thr | sim(:)>thr)));

HIT = @(sim,obs)D1H1(sim,obs)./(D0H1(sim,obs)+D1H1(sim,obs));
FALSE_ALARM = @(sim,obs)D1H0(sim,obs)./(D1H0(sim,obs)+D1H1(sim,obs));
CRITICAL = @(sim,obs)D1H1(sim,obs)./(D1H1(sim,obs)+D0H1(sim,obs)+D1H0(sim,obs));
ERROR_BIAS = @(sim,obs)D1H0(sim,obs)./D0H1(sim,obs);

STATS.rmse = NaN(size(FLO_test.floodmaps,1),1);

mapsSim = FLO_pred.floodmaps;
mapsObs = FLO_test.floodmaps;
STATS.eventTi = [];


mapSim = getmapSim(mapsSim,FLO_pred.eventNo,DIS_analog,testConfig);


for evi = 1:size(mapsSim,1)
    
    
    STATS.rmse(evi) = RMSE(aggregateImage(mapSim(evi,:,:),aggUnit,'max'),...
        aggregateImage(mapsObs(evi,:,:),aggUnit,'max'));
    STATS.mape(evi) = MAPE(aggregateImage(mapSim(evi,:,:),aggUnit,'max'),...
        aggregateImage(mapsObs(evi,:,:),aggUnit,'max'));
    
    STATS.eventTi(evi) = evi - find(FLO_test.eventNo == FLO_test.eventNo(evi),1);
    STATS.eventTT(evi) = find(FLO_test.eventNo == FLO_test.eventNo(evi),1,'last')-...
        find(FLO_test.eventNo == FLO_test.eventNo(evi),1)+1;
    STATS.hit(evi) = HIT(aggregateImage(mapSim(evi,:,:),aggUnit,'max'),...
        aggregateImage(mapsObs(evi,:,:),aggUnit,'max'));
    STATS.fa(evi) = FALSE_ALARM(aggregateImage(mapSim(evi,:,:),aggUnit,'max'),...
        aggregateImage(mapsObs(evi,:,:),aggUnit,'max'));
    STATS.cci(evi) = CRITICAL(aggregateImage(mapSim(evi,:,:),aggUnit,'max'),...
        aggregateImage(mapsObs(evi,:,:),aggUnit,'max'));
    STATS.eb(evi) = ERROR_BIAS(aggregateImage(mapSim(evi,:,:),aggUnit,'max'),...
        aggregateImage(mapsObs(evi,:,:),aggUnit,'max'));
    
end

STATS.eventNo = FLO_test.eventNo;
STATS.eventNo_pred = FLO_pred.eventNo;

end



function mapSim = getmapSim(mapsSim,eventNoSim,DIS_analog,testConfig)

if testConfig.uniqueEvent == false && testConfig.IDW == false
    mapSim = squeeze(nanmedian(mapsSim,4));
    if testConfig.includePrevious == true
        mapSim = cumsum(mapSim,1)./reshape(1:size(mapSim,1),[size(mapSim,1),1,1]);
    end
    
    return
end


if testConfig.uniqueEvent == false
    
    if testConfig.IDW == true
        
        mapSim = NaN(size(mapsSim,1:3));
        
        for evti = 1:size(mapsSim,1)
            wei = getIDWweight(DIS_analog(evti,:),1);
            oriSim = squeeze(mapsSim(evti,:,:,:));
            mapSim(evti,:,:) = squeeze(nansum(oriSim.*reshape(wei,[1,1,length(wei)]),3));
        end
        
    end
    
else
    
    mapSim = NaN(size(mapsSim,1:3));
    
    for evti = 1:size(mapsSim,1)
        
        oriSimEventNo = eventNoSim(evti,:);
        
        oriSim = squeeze(mapsSim(evti,:,:,:));
        [~,uniI,~] = unique(oriSimEventNo);
        mapSim(evti,:,:) = squeeze(nanmedian(squeeze(oriSim(:,:,uniI)),3));
        % mapSim(evti,:,:) = squeeze(nansum(oriSim.*reshape(wei,[1,1,length(wei)]),3));
    end
    
end

if testConfig.includePrevious == true
    mapSim = cumsum(mapSim,1)./reshape(1:size(mapSim,1),[size(mapSim,1),1,1]);
end

end


