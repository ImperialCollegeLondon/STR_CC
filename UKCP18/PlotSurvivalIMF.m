%% Configuration
Name = {'London'};%{'London','SWestuk','Westuk','Scotland'};

%%
CPM = load(sprintf('rawIMF_Area_%s_4Season.mat',Name{thisstation}),'rawIMF');
RAD = load(sprintf('Obs_rawIMF_Area_%s_4Season.mat',Name{thisstation}),'rawIMF');

warning off

seasonNo = 3;
aggNo = 20;
enNo = 1;
colm = pink(8);
colm(1:3,:) = [];
colObs = gray(5);
% setFigureProperty;
sptag = 1;
for seasonNo = 1:4
for aggNo = [1 5 10 24]
    subplot(4,4,sptag)
    aggArea = CPM.rawIMF{seasonNo}(enNo).areaAgg(aggNo);
    for enNo = 1:12
        [f,x,flo,fup] = ecdf(CPM.rawIMF{seasonNo}(enNo).AreaIMFs{aggNo});
        h = plot(x,1-f,'color',colm(1,:),'linewidth',1)
        hold on;
    end
    [f,x,flo,fup] = ecdf(RAD.rawIMF{seasonNo}.AreaIMFs{aggNo});
    hobs = plot(x,1-f,'color',colObs(2,:),'linewidth',2);
    % set(gca,'XScale','log')
    
    set(gca,'YScale','log')
    XL = xlim;
    xlim([2,XL(2)])
    grid minor
    YL = ylim;
    ylim([10^-5,YL(2)])
    text(2.5,YL(2)*0.50,{sprintf('%4.2f Km^2',aggArea);getSeasonName(seasonNo)});
    legend([h,hobs],{'R1-UKCP','R1-Radar'});
    xlabel('Areal Mean Rainrate mm/h')
    ylabel('Exceedence probability');
    sptag = sptag+1;
end
end
XYWH = [50,-50,500,500];
set(gcf,'units','points','position',XYWH);

