clear;clc

addpath(genpath('C:\Users\Yuting Chen\Dropbox (Personal)\MatLab_Func'));
addpath(genpath(cd));


REGIONS = REGIONS_info();
dataSP = 'H:\DATA_CLIMATE\UKCP18\';
region = REGIONS.Westuk;%London;%Scotland;%SWestuk;%
seasonNo = 2;


%%
figure;
XYWH = [50,-50,500,360];
set(gcf,'units','points','position',XYWH);
ha = tight_subplot(2,3,[.01 .01],[.1 .01],[.1 .01]);

h = [];
colm = pink(15);
linewidth = 1.5;
linewidthObs = 3;
% legend(h,)
OBS = load(sprintf('%sCellProp_NIMROD_Season%01d_%s.mat',dataSP,seasonNo,region.Name),...
    'CELLA','THRE');
for i = 1:length(OBS.THRE)
    thisplot = i;
    axes(ha(thisplot));
    aa = OBS.CELLA{i}{1};
    [f,x,flo,fup] = ecdf(aa);
    h(thisplot) = plot(x*2.2^2,1-f,':','color',colm(i,:),'Linewidth',linewidth);
    hold on;
end

CPM = load(sprintf('%sCellProp_UKCP_Season%01d_%s.mat',dataSP,seasonNo,region.Name),...
    'CELLA','THRE');
for i = 1:length(CPM.THRE)
    thisplot = i;
    axes(ha(thisplot));
    for enNo = 1:12
        [f,x,flo,fup] = ecdf(CPM.CELLA{i}{enNo});
        hsim(thisplot) = plot(x*2.2^2,1-f,'-','color',colm(i,:),'Linewidth',linewidthObs);
        hold on;
    end
    % set(gca,'XScale','log');
    xlabel('Cell Area Km^2')
    ylabel('1-F(x)');
    xlim([0^2,110*2.2^2])
    grid minor
    set(gca,'XGrid','on');
    text(6,0.05,sprintf('%dmm/h',CPM.THRE(i)),'fontsize',12);
    legend(gca,[h(thisplot),hsim(thisplot)],...
                {'2007-2017 Radar','UKCP 12Ensembles'},'Location','NorthEast');
end
format_xylabel(ha,2,3)


