% THIS SCRIPT IS TO INVESTIGATE HOW TO INCLUDE THE PREVIOUS IMAGES

function aux_modelSetting(FOREOUTPUT,STATS)

% load('K:\DATA_FCS\CrossVal\v1_PCA_PCs05\STATS_v1_ev004.mat')

testConfig = struct();
testConfig.IDW = false;% false is better (at least for event 90)
testConfig.uniqueEvent = false;
testConfig.includePrevious = false;

tic
[STATS] = evaluateFloodPrediction(FOREOUTPUT.FLO_test,FOREOUTPUT.FLO_pred,...
    FOREOUTPUT.E,FOREOUTPUT.N,...
    FOREOUTPUT.analogDistL1,testConfig);
toc

close all
ha = tight_subplot(1,4,[.05 .05],[.1 .1],[.15 .05]);
axes(ha(1))
imagesc(FOREOUTPUT.analogDistL1);shading flat
cptcmap('illumination', 'mapping','scaled','flip',false);%,'ncol',20);
xlabel('Ensemble No.')
ylabel('Time to storm origin')
caxis([0,2])
colorbar('SouthOutside')
set(gca,'YDir','normal')
title('Seuclidean Distance Layer 1')

axes(ha(2))
imagesc(FOREOUTPUT.analogDistL2);shading flat
cptcmap('illumination', 'mapping','scaled','flip',false);%,'ncol',20);
xlabel('Ensemble No.')
% ylabel('Time to storm origin')
caxis([0,10])
colorbar('SouthOutside')
set(gca,'YDir','normal')
title('Seuclidean Distance Layer 2')

axes(ha(3))
stats = [STATS.hit;1-STATS.fa;STATS.cci;STATS.eb];
imagesc(stats');shading flat
cptcmap('diff_4_bias', 'mapping','scaled','flip',false,'ncol',25);
caxis([-3.5,5.5])
ax = gca;
ax.XTick = 1:4;
ax.XTickLabel = {'HR','FA','CSI','EB'};
% ylabel('Time to storm origin')
xlabel('Evaluation Metrics')
set(gca,'YDir','normal')
c = colorbar('SouthOutside');
c.Ticks = [-3.5,1,5.5];
c.TickLabels = {'Bad','Perfect','Bad'};     
title('Evaluation Metrics')

axes(ha(4))
imagesc(FOREOUTPUT.analogPredEvNo);shading flat
cptcmap('GMT_gray', 'mapping','scaled','flip',false,'ncol',157);
caxis([1,157])
title('Identified Event')
set(gca,'YDir','normal')
colorbar('SouthOutside')
% ylabel('Time to storm origin')
xlabel('Ens No.')

end
