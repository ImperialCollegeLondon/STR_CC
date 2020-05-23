
close all;
evNo = '013';

figure(1)
setFigureProperty()
load(['K:\DATA_FCS\CrossVal\v1_PCA_PCs05\STATS_v1_ev',evNo,'.mat'])
aux_modelSetting(FOREOUTPUT,STATS)
XYWH = [0,0,1400,400];
set(gcf,'units','points','position',XYWH);

figure(2)
load(['K:\DATA_FCS\CrossVal\v1_PCA_PCs05_RainRAD\STATS_v1_ev',evNo,'.mat'])
aux_modelSetting(FOREOUTPUT,STATS)
XYWH = [0,350,1400,400];
set(gcf,'units','points','position',XYWH);








