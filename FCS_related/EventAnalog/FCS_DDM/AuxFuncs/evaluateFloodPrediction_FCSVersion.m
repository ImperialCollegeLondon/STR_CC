function [STATS] = evaluateFloodPrediction_FCSVersion(FLO_test,FLO_pred,E,N,DIS_analog,testConfig)
% this function is to compute several evaluation metrics as we agreed on
% FCS meeting.

STATS = struct;
unit = '3km';
aggUnit = 3;
thr = 0;


TP = @(sim,obs)1e-6+sum(sim(:)>thr & obs(:)>thr);
FN = @(sim,obs)1e-6+sum(sim(:)<=thr & obs(:)>thr);
FP = @(sim,obs)1e-6+sum(sim(:)>thr & obs(:)<=thr);
TN = @(sim,obs)1e-6+sum(sim(:)<=thr & obs(:)<=thr);


RMSE = @(sim,obs)sqrt(nanmean((squeeze(sim(:))-squeeze(obs(:))).^2));
MAPE = @(sim,obs)nanmean(abs(sim(obs(:)>thr | sim(:)>thr) - ...
    obs(obs(:)>thr | sim(:)>thr)./obs(obs(:)>thr | sim(:)>thr)));

TPR = @(sim,obs)TP(sim,obs)./(TP(sim,obs)+FN(sim,obs));
TNR = @(sim,obs)TN(sim,obs)./(FP(sim,obs)+TN(sim,obs));
PPR = @(sim,obs)TP(sim,obs)./(TP(sim,obs)+FP(sim,obs));
NPR = @(sim,obs)TN(sim,obs)./(TN(sim,obs)+FN(sim,obs));
ACC = @(sim,obs)(TN(sim,obs)+TP(sim,obs))./...
    (TP(sim,obs)+TN(sim,obs)+FN(sim,obs)+FP(sim,obs));

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
    STATS.tpr(evi) = TPR(aggregateImage(mapSim(evi,:,:),aggUnit,'max'),...
        aggregateImage(mapsObs(evi,:,:),aggUnit,'max'));
    STATS.tnr(evi) = TNR(aggregateImage(mapSim(evi,:,:),aggUnit,'max'),...
        aggregateImage(mapsObs(evi,:,:),aggUnit,'max'));
    STATS.ppr(evi) = PPR(aggregateImage(mapSim(evi,:,:),aggUnit,'max'),...
        aggregateImage(mapsObs(evi,:,:),aggUnit,'max'));
    STATS.npr(evi) = NPR(aggregateImage(mapSim(evi,:,:),aggUnit,'max'),...
        aggregateImage(mapsObs(evi,:,:),aggUnit,'max'));
    STATS.acc(evi) = ACC(aggregateImage(mapSim(evi,:,:),aggUnit,'max'),...
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


