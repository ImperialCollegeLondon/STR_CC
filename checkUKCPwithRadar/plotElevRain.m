% plot Elevation - Rain





%% get CROSS_CPM_RADAR_GEAR

for mon = 1:12
    
    figure;
    ha = tight_subplot(4,1,[.02 .02],[.1 .05],[.1 .1]);
    
    vq_cpm = [];
    for EnsNo = 1:12
        F_cpm = scatteredInterpolant(E(:)*1000,N(:)*1000,CPM{mon,EnsNo}(:),'nearest');
        vq_cpm{EnsNo} = F_cpm(pCross(:,1),pCross(:,2));
    end
    
    F_rad = scatteredInterpolant(E(:)*1000,N(:)*1000,squeeze(RAD(mon,:)'),'nearest');
    vq_rad = F_rad(pCross(:,1),pCross(:,2));
    
    
    F_gear = scatteredInterpolant(X_coor(:)*1000,Y_coor(:)*1000,squeeze(GEAR(mon,:)'),'nearest');
    vq_gear = F_gear(pCross(:,1),pCross(:,2));
    
    axes(ha(1))
    plot(vq)
    formatItCrossElev('Terrain50',mon,'Topo');
    
    axes(ha(2))
    area(vq_cpm{1},'FaceColor','r')
    formatItCrossElev('2.2km (1980-2000)',mon,'Ens');
    for EnsNo = 2:12
        area(vq_cpm{EnsNo},'FaceColor','r')
        formatItCrossElev([],mon,'Ens')
    end
    legend('Ensembles')
    grid minor
    
    axes(ha(3))
    area(vq_rad,'FaceColor','b')
    formatItCrossElev('Radar (2007-2018)',mon)
    
    axes(ha(4))
    plot(vq_gear/24/eomday(2000,mon),'-');hold on
    vq_gear = smooth(vq_gear,3);
    area(vq_gear/24/eomday(2000,mon),'FaceColor','b')
    formatItCrossElev('GEAR (1980-2000)',mon)
    legend('1km grid','~2.2km')
    xlabel('Cross-section')
    
    set(ha(1:3),'XTickLabel',[]);
    
    
    filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP';
    filename = [filePath,filesep,'Topo_Rainfall',sprintf('%02d',mon)];
    savePlot(filename,'XYWH',[150,0,400,700],'NeedReply','N');
    
    close all
    
end




function formatItCrossElev(textStr,mon,variable)
arguments
    textStr (1,:) char
    mon (1,1) double
    variable (1,:) char = 'Rain'
end
hold on;
xlim([1,68])
if strcmpi(variable,'Rain')
    alpha(0.4)
    YLIM = ylim;
    ylim([0,ceil(YLIM(2)*1.25)/1.25])
end

if strcmpi(variable,'Ens')
    alpha(0.2)
    YLIM = ylim;
    ylim([0,ceil(YLIM(2)*1.25)/1.25])
end


grid minor
YLIM = ylim;
if ~isempty(textStr)
    text(1.5,YLIM(2),sprintf('%s-%s',textStr,getMonthName(mon,true)),'VerticalAlignment','Top')
end
end