

fprintf('-----------------------\n');
season = 'JJA';
option = 'm2Hourly';%'thre50Hourly';% 'm2Daily';%
% CONFIGURATION
obscol = 'k';
simcol = 'r';
p = 75;
minVal = -Inf;%-Inf;%Inf;%-0.5;
% PLOT
figure;
setFigureProperty('Paper');
hold on
filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP';

switch(option)
    case 'm2Hourly'
        load(['D:\UKCP18\Radar_4d4_JJA_pJJAp2_RipletK2.mat'])
        repInd = [];
        filename = [filePath,filesep,'BesagL_',season,'_p',season,'p2'];% figsave name
    case 'thre50Hourly'
        load('D:\UKCP18\Radar_4d4_JJA_Thre50_RipletK2.mat')
        repInd = 6:10;% somehow I compute several radius for twice???
        filename = [filePath,filesep,'BesagL_',season,'_thre50'];% figsave name
    case 'm2Daily'
        load('D:/UKCP18/Radar_4d4_JJA_pJJADaily_p2_RipletK2.mat');
        repInd = [];
        filename = [filePath,filesep,'BesagL_',season,'_m2Daily'];% figsave name
end

x = RadVec*4.4;
y = sqrt(AK/pi)-RadVec;
x(repInd) = [];RadVec(repInd) = []; y(:,repInd) = [];
y = y(nanmin(y,[],2)>=minVal,:);
hobs = plot(x,nanmean(y,1),'-','color',obscol,'linewidth',3);
% plot(x,prctile(y,[p],1),'--','color',obscol,'linewidth',2);
% plot(x,prctile(y,[100-p],1),'--','color',obscol,'linewidth',2);
fprintf('Radar Observation:%03f/year\n',size(AK,1)/12);

ENSEMBLENO=getEnsNos();
[AK] = deal([]);
for eni = 1:length(ENSEMBLENO)
    try
        switch(option)
            case 'm2Hourly'
                C = load(['D:\UKCP18\CPM_4d4_Ensno',ENSEMBLENO{eni},'_',season,...
                    '_p',season,'p2_RipletK2.mat']);%
            case 'thre50Hourly'
                C = load(['D:\UKCP18\CPM_4d4_Ensno',ENSEMBLENO{eni},'_',season,...
                    '_Thre50_RipletK2.mat'])';
            case 'm2Daily'
                C = load(['D:\UKCP18\CPM_4d4_Ensno',ENSEMBLENO{eni},'_',season,...
                    '_p',season,'Daily_p2_RipletK2.mat']);%
                
        end
        C.AK(end:2000,:) = NaN;
        
        AK = cat(3,AK,reshape(C.AK,[size(C.AK),1]));
        RadVec = C.RadVec;
        y = sqrt(AK/pi)-RadVec;
    catch me
        break
    end
end

RadVec(repInd) = []; y(:,repInd,:) = [];
for ens = 1:size(y,3)
    y(nanmin(squeeze(y(:,:,ens)),[],2)<=minVal,:,ens) = NaN;
end
sm = @(x)reshape(smooth(x,1),[1,numel(x)]);
y = squeeze(nanmean(y,1));
hsim = plot(x,nanmedian(y,2),'-','color',simcol,'linewidth',2);
hsimran = fill([x,flip(x)],[sm(nanmin(y,[],2)),sm(flip(nanmax(y,[],2)))],'r',...
    'LineStyle','none');
alpha(0.3)
% plot(x,prctile(y,[p],1),'--','color',simcol,'linewidth',1);
% plot(x,prctile(y,[100-p],1),'--','color',simcol,'linewidth',1);
fprintf('CPM EnsNo%s:%03f/year\n',ENSEMBLENO{eni},size(C.AK,1)/20);

hold off

% FORMAT
ylabel('Besag L');
xlabel('r(L)(km)');
% grid minor
box on
legend([hobs,hsim,hsimran],{'RADAR','CPM2.2km','Ensemble Range'})
set(gca,'linewidth',2)
ylim([0,20])

fprintf('-----------------------\n');

% SAVE
savePlot(filename,'units','centimeters','XYWH',[5,0,8,6],'needreply','N');

