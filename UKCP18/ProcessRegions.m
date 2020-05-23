% Process regional Pr
% London Region;

%%
clear;clc

REGIONS = REGIONS_info();
dataSP = 'H:\DATA_CLIMATE\UKCP18\';

Period = '2060-2080';%'1980-2000';%
region = REGIONS.London;%SWestuk;%Westuk;%Scotland;
ProcessIt(region,dataSP,Period);
fprintf('%s Finished.\n',region.Name);

region = REGIONS.SWestuk;%London;%Westuk;%Scotland;
ProcessIt(region,dataSP,Period);
fprintf('%s Finished.\n',region.Name);

region = REGIONS.Westuk;%London;%SWestuk;%Scotland;
ProcessIt(region,dataSP,Period);
fprintf('%s Finished.\n',region.Name);

region = REGIONS.Scotland;%London;%SWestuk;%Westuk;%
ProcessIt(region,dataSP,Period);
fprintf('%s Finished.\n',region.Name);
tag = 1;

%%


function ProcessIt(region,dataSP,Period)

%% FIGURE1. Save data to plot WetArea(Thres) from CPM and Radar respectively

% for cpm

[dataplot] = cell(1,4);

for thissp = 1:4
    
    season = thissp;
    
    [E,N,RainEnsembles,scaleF,region] = loadRainEnsembles(region,getMons(season),Period);
    [Thres,y1,y2] = compute_WetArea(E,N,RainEnsembles,scaleF,region);
    dataplot{thissp} = struct('Thres',Thres,'y1',y1,'y2',y2); %#ok<SAGROW>
    
    fprintf('------thissp_%1d------\n',thissp)
    
end

save(sprintf('%sWAR_Area_%s_%s_4Season.mat',dataSP,region.Name,Period),'dataplot');
clear RainEnsembles


% for radar
[dataplot] = cell(1,4);

for thissp = 1:4
    
    season = thissp;

    [E,N,RainNimrod,scaleF,region] = loadRainNimrod(region,getMons(season));
    [Thres,y1,y2] = compute_WetArea(E,N,{RainNimrod},scaleF,region);
    dataplot{thissp} = struct('Thres',Thres,'y1',y1,'y2',y2); %#ok<SAGROW>
    
    fprintf('------thissp_%1d------\n',thissp)
    
end
save(sprintf('%sObs_WAR_Area_%s_4Season.mat',dataSP,region.Name),'dataplot');


%% FIGURE2. Save data to plot Intensity(Area)


% for CPM
dataplot = cell(1,4);
rawIMF = [];
for thissp = 1:4
    seasonInd = thissp;
    [E,N,RainEnsembles,scaleF,region] = loadRainEnsembles(region,getMons(seasonInd),Period);
    [areaAgg,y1,y2,AreaIMFs] = compute_IMF(E,N,RainEnsembles,scaleF,region);
    clear RainEnsembles;
    dataplot{thissp} = struct('areaAgg',areaAgg,'y1',y1,'y2',y2);
    rawIMF{thissp} = struct('AreaIMFs',AreaIMFs,'areaAgg',areaAgg);
    fprintf('------thissp_%1d------\n',thissp) 
end
save(sprintf('%sIMF_Area_%s_%s_4Season.mat',dataSP,region.Name,Period),'dataplot');
save(sprintf('%srawIMF_Area_%s_%s_4Season.mat',dataSP,region.Name,Period),'rawIMF')

% for radar
[dataplot] = cell(1,4);
for thissp = 1:4
    season = thissp;
    [E,N,RainNimrod,scaleF,region] = loadRainNimrod(region,getMons(season));
    [areaAgg,y1,y2,AreaIMFs] = compute_IMF(E,N,{RainNimrod},scaleF,region);
    dataplot{thissp} = struct('areaAgg',areaAgg,'y1',y1,'y2',y2);
    rawIMF{thissp} = struct('AreaIMFs',AreaIMFs,'areaAgg',areaAgg);
    fprintf('------thissp_%1d------\n',thissp)
end
save(sprintf('%sObs_IMF_Area_%s_4Season.mat',dataSP,region.Name),'dataplot');
save(sprintf('%sObs_rawIMF_Area_%s_4Season.mat',dataSP,region.Name),'rawIMF')


%% FIGURE3. Save data to plot Intensity-Duration (Area-based)

% for cpm

[dataplot] = cell(1,4);

for thissp = 1:4
    
    season = thissp;
    
    [E,N,RainEnsembles,scaleF,region] = loadRainEnsembles(region,getMons(season),Period);
    [~,~,WAR] = compute_WetArea(E,N,RainEnsembles,scaleF,region,0.1);
    for ensNo = 1:length(WAR)
        WAR0 = WAR{ensNo}{1};
        IMF = computeIMF(RainEnsembles{ensNo},scaleF,region.gridArea);
        PIMF = nanmax(reshape(double(RainEnsembles{ensNo})/scaleF,[],size(RainEnsembles{ensNo},3)),[],1);
        [lengthi,stormsi] = seperateRadarEvents( WAR0, 0.02, 2, 60, 60);
        EvNo = [zeros(1,stormsi(1)-1),cell2mat(arrayfun(@(si,nsi,li,evi)...
            [evi*ones(1,li),0*ones(1,nsi-si-li)],stormsi,...
            [stormsi(2:end),length(IMF)+1],lengthi,1:length(stormsi),'UniformOutput', false))];
        T = table(IMF(EvNo>0)',PIMF(EvNo>0)',EvNo(EvNo>0)',WAR0(EvNo>0)',...
            'VariableNames',{'IMF','PIMF','EvNo','WAR'});
        eIMF = grpstats(T,'EvNo','mean').mean_IMF;
        eWAR = grpstats(T,'EvNo','median').median_WAR;
        pIMF = grpstats(T,'EvNo','max').max_PIMF;
        pWAR = grpstats(T,'EvNo','max').max_WAR;
        dataplot{thissp}{ensNo} = struct('IMF',eIMF,'pIMF',pIMF,'Dur',lengthi',...
            'WAR',eWAR,'pWAR',pWAR,'stormsi',stormsi'); %#ok<SAGROW>
        
    end
    % dataplot{thissp} = struct('Thres',Thres,'y1',y1,'y2',y2); %#ok<SAGROW>
    
    fprintf('------thissp_%1d------\n',thissp)
    
end

save(sprintf('%sDUR_Area_%s_%s_4Season.mat',dataSP,region.Name,Period),'dataplot');


%  for radar
[dataplot] = cell(1,4);
for thissp = 1:4
    season = thissp;
    [E,N,RainNimrod,scaleF,region] = loadRainNimrod(region,getMons(season));
    [~,~,WAR] = compute_WetArea(E,N,{RainNimrod},scaleF,region,0.1);
    WAR = WAR{1}{1};
    IMF = computeIMF(RainNimrod,scaleF,region.gridArea);
    PIMF = nanmax(reshape(double(RainNimrod)/scaleF,[],size(RainNimrod,3)),[],1);
    % saveFigFp = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP';
    % hold on
    % bar((1:length(WAR)),WAR);plot((1:length(WAR)),WAR,'k');
    % plot(stormsi,(stormsi>0)-1,'r.','markersize',50);
    % plot((stormsi+lengthi-1),(stormsi>0)-1,'k.','markersize',50);
    % hold off;box on;xlabel('Time (Hour)');ylabel('WAR(>0.02)');title('Event Seperation')
    % set(gca,'linewidth',2)
    % figName = 'TemporalAnalysis\EventSepration';
    % % savePlot([saveFigFp,filesep,figName],'XYWH',[50,50,500,220],'needreply','N','onlyPng',true);
    [lengthi,stormsi] = seperateRadarEvents( WAR, 0.02, 2, 60, 60);
    EvNo = [zeros(1,stormsi(1)-1),cell2mat(arrayfun(@(si,nsi,li,evi)...
        [evi*ones(1,li),0*ones(1,nsi-si-li)],stormsi,...
        [stormsi(2:end),length(IMF)+1],lengthi,1:length(stormsi),'UniformOutput', false))];
    T = table(IMF(EvNo>0)',PIMF(EvNo>0)',EvNo(EvNo>0)',WAR(EvNo>0)',...
            'VariableNames',{'IMF','PIMF','EvNo','WAR'});
    eIMF = grpstats(T,'EvNo','mean').mean_IMF;
    eWAR = grpstats(T,'EvNo','median').median_WAR;
    pIMF = grpstats(T,'EvNo','max').max_PIMF;
    pWAR = grpstats(T,'EvNo','max').max_WAR;
    dataplot{thissp} = struct('IMF',eIMF,'pIMF',pIMF,'Dur',lengthi',...
        'WAR',eWAR,'pWAR',pWAR,'stormsi',stormsi'); %#ok<SAGROW>
end

save(sprintf('%sObs_DUR_Area_%s_4Season.mat',dataSP,region.Name),'dataplot');


end

% %%
% 
% figure;
% for enNo = 1:12
%     subplot(3,4,enNo);
%     monMean = 365*24*getSpaceMean(RainEnsembles{enNo},scaleF);
%     contourf(E,N,monMean,[600:200:3000]);
%     col = winter(20);
%     col = parula(20);
%     colormap(col(end-6:-1:4,:))
%     title(sprintf('Ensemble No.%02d',getEnsembleNo(enNo)));
%     axis off
%     box on
% end
% colorbar
% 

%% Power Spectrum for one sites
% figure;
% [y1,y2] = deal([]);
% 
% figure;
% Pf = [];
% itag = 1;
% for enNo = 1:length(RainEnsembles)
%     
%     R3d = double(RainEnsembles{enNo})/scaleF;
%     y1{enNo} = squeeze(R3d(:,:,1));
%     y2{enNo} = squeeze(R3d(25,25,:));
%     
%     RAPs_tem = squeeze(R3d(25,25,:));
%     RAPs_tem = RAPs_tem(~isnan(RAPs_tem));%RAPs_tem((RAPs_tem>0));
%     % Plot RAPS
%     pl = 0;
%     col = [0.5 0.5 0.5];
%     for i = 1:floor(length(RAPs_tem)/1000)
%         try
%             [f1,Pf{itag},xCell,yCell] = raPsd2d(RAPs_tem(1000*(i-1)+(1:1000)),1,pl);
%             plot_RAPS(f1,Pf{itag},xCell,yCell,col)
%             itag = itag + 1;
%         catch
%             1;
%         end
%     end
%     
% end
% medianPf = nanmean(cell2mat(Pf'),1);
% colMedian = 'r';
% plot_RAPS(f1,medianPf,xCell,yCell,colMedian);
% ylim([1e-8,1]);
% RAPS = struct('E',E(25,25),'N',N(25,25),'Ensembles','All12');
% save(sprintf('%sRAPS_%s_mon%-2d.mat',dataSP,region.Name,mon),'f1','Pf','RAPS');

%% several necessary function
% load group



function imf = computeIMF(m3d,scaleF,areaMat)

if isvector(m3d)
    imf = double(m3d)/scaleF;
else
    m3d = squeeze(m3d);
    imf = nansum(reshape(double(m3d)/scaleF.*areaMat,[],size(m3d,3)),1)...
        ./nansum(areaMat(:));
    % here, 'convn' is also OK!
end
end

% get group
function ensemNo = getEnsembleNo(ind)

ENSEMBLENO=[1,4,5,6,7,8,9,10,11,12,13,15];
ensemNo = ENSEMBLENO(ind);

end


function mR = getSpaceMean(Mat,scaleF)

mR = squeeze(nanmean(Mat,3)/scaleF);

end

function [E1,N1] = getEN(region)
filePath = ['J:/UkCp18/','01','/pr_rcp85_land-cpm_uk_2.2km_'];
ncFileName = [filePath,'01','_1hr_19910801-19910830.nc'];
LAT=ncread(ncFileName,'latitude');
LON=ncread(ncFileName,'longitude');
[E, N] = ll2os(LAT, LON);
E = E/1000;
N = N/1000;
[region.i,region.j] = getRegionIJ(E,N,region.minE,region.minN);

E1 = E(region.i:region.i+region.dimE-1, region.j);
N1 = N(region.i, region.j:region.j+region.dimN-1);

[E1,N1] = meshgrid(E1,N1);
E1 = E1';
N1 = N1';

end

% compute statistcis group
function [Thres,meanWAR,WAR] = compute_WetArea(E,N,RainEnsembles,scaleF,region,Thres)
arguments
    E double
    N double
    RainEnsembles cell
    scaleF double
    region struct
    Thres (1,:) double = [1:4,5:2:25]
end
% figure;
colm = pink(12);
[meanWAR,WAR] = deal([]);

for enNo = 1:length(RainEnsembles)
    
    [allWAR] = arrayfun(@(thre)(computeWAR(RainEnsembles{enNo},thre,scaleF,region)),...
        Thres, 'UniformOutput', false);
    meanWAR(enNo,:) = cell2mat(cellfun(@(x)nanmean(x),allWAR,'UniformOutput', false));
    WAR{enNo} = allWAR;% cell2mat(cellfun(@(x)x,WAR,'UniformOutput', false));
    % war = computeWAR(RainEnsembles{enNo},thre,scaleF);
    
end

end
% compute SUB-GROUP
function [war] = computeWAR(m3d,thre,scaleF,region)
% considering the grid area is not perfectly 2.2km, so here we also
% included the grid area size for each grid.

war = nanmean(reshape((m3d>thre*scaleF).*region.gridArea,[],size(m3d,3)),1)...
    ./nanmean(double(region.gridArea(:)));

end

function [areaAgg,y1,y2,AreaIMFs] = compute_IMF(E,N,RainEnsembles,scaleF,region)
[y1,y2,AreaIMFs] = deal([]);

for enNo = 1:length(RainEnsembles)
    
    lenAgg = 1:24;%(1:50).^2;
    
    % imf = computeIMF(m3d,thre,scaleF,areaMat);
    
    ranAgg = @(lenAgg)25-lenAgg:25+lenAgg;
    areaSize = @(lenAgg)region.gridArea(ranAgg(lenAgg),ranAgg(lenAgg));
    eachIMF = @(lenAgg)computeIMF(RainEnsembles{enNo}...
        (ranAgg(lenAgg),ranAgg(lenAgg),:),scaleF,areaSize(lenAgg));
    
    AreaIMF = arrayfun(eachIMF,...
        lenAgg, 'UniformOutput', false);
    AreaIMFs{enNo} = AreaIMF;
    
    y1(enNo,:) = cell2mat(cellfun(@(x)nanmean(x),AreaIMF,'UniformOutput', false));
    y2(enNo,:) = cell2mat(cellfun(@(x)nanmean(x(x>0)),AreaIMF,'UniformOutput', false));
    
    % war = computeWAR(RainEnsembles{enNo},thre,scaleF);
    
end
areaAgg = cell2mat(arrayfun(@(lenAgg)nansum(nansum(region.gridArea...
    (ranAgg(lenAgg),ranAgg(lenAgg)))),...
        lenAgg, 'UniformOutput', false));
end

function plot_RAPS(f1,Pf,xCell,yCell,col)

fontSize = 14;
% setFigureProperty;

loglog(1./f1,Pf,'-','LineWidth',0.5,'color',col)
set(gcf,'color','white')
set(gca,'FontSize',fontSize,'FontWeight','bold',...
    'XGrid','on','YAxisLocation','right','XDir','reverse');
% set(gca,'YTickLabel',yCell,'YMinorTick','off',...
%     'XTickLabel',xCell)
xlabel('Wavelength (h)','FontSize',fontSize,'FontWeight','Bold');
ylabel('Power','FontSize',fontSize,'FontWeight','Bold');
title('Radially averaged power spectrum','FontSize',fontSize,'FontWeight','Bold')
hold on;

end











