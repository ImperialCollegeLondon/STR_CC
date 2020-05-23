clear;clc
close all

durVec = [1 3 6 12 24];

for regionName = {'EUK','WAL','SCO'}%{'WAL'}%
    
    regionName = regionName{1};
    
    dataSourcePath = 'D:\UKCP18\fixedAreaARF';%'D:\UKCP18\stormcentredARF';%
    close all
    figure;
    setFigureProperty('Subplot2')
    
    dti = 1;% 1:length(durVec)
    dur = durVec(dti);
    colpas = getColorMat(dti,'pas');
    colfut = getColorMat(dti,'fut');
    
    if ~isscalar(dti)
        error(1);
    end
    ensInd  = 0;
    [allRatio,arf_pas,arf_fut] = deal([]);
    for ensNo = getEnsNos()
        [Dat_past,Dat_future,arf] = getData(dataSourcePath,regionName,ensNo,dur);
        % [Dat_past,Dat_future,arf] = getData(dataSourcePath,regionName,getEnsNos(),dur);
        
        % subplot 121
        hold on;
        x = [1,arf.area]*2.2*2.2;
        y = [1,nanmean(Dat_past.arf,1)];
        %{
        hpre = plot(x,y,...
            '--','color',[colpas,0.3],'markerfacecolor',[colpas],'markersize',5);
        y = [1,nanmean(Dat_future.arf,1)];
        hfut = plot(x,y,...
            '--','color',[colfut,0.3],'markerfacecolor',[colfut],'markersize',10);
        hold off
        %}
        %{
        subplot 122
        hold on;
        % [hpre,hfut] = deal(handle(0));
        col = getColorMat(dti,'fut');
        x = [arf.area]*2.2*2.2;
        y = nanmean(Dat_future.arf,1)./...
            nanmean(Dat_past.arf,1);
        hpre = plot(x,y,...
            '-','color',[col,0.3],'markerfacecolor',[col]);
        hold off
        %}
        ensInd = ensInd+1;
        allRatio(ensInd,:) = nanmean(Dat_future.arf,1)./...
            nanmean(Dat_past.arf,1);
        arf_pas(ensInd,:) = [1,nanmean(Dat_past.arf,1)];
        arf_fut(ensInd,:) = [1,nanmean(Dat_future.arf,1)];
        arf_x = [1,arf.area]*2.2*2.2;
        allRatio_x = [arf.area]*2.2*2.2;
        
        drawnow
        
    end
    %%
    
    % subplot 121
    hold on;
    
    hpre = plot(arf_x,nanmedian(arf_pas,1),...
        '-','color',[colpas,0.8],'markerfacecolor',[colpas],'markersize',5);
    fill([arf_x,flip(arf_x)],[nanmax(arf_pas,[],1),flip(nanmin(arf_pas,[],1))],...
        [colpas],'LineStyle','none');
    hfut = plot(arf_x,nanmedian(arf_fut,1),...
        '-','color',[colfut,0.8],'markerfacecolor',[colfut],'markersize',10);
    fill([arf_x,flip(arf_x)],[nanmax(arf_fut,[],1),flip(nanmin(arf_fut,[],1))],...
        [colfut],'LineStyle','none');
    alpha(0.3)
    % formatSubP1(hpre,hfut,durVec([1]))
    hold off
    
    yyaxis right
    hold on;
    col = getColorMat(dti,'fut');
    hrat = plot(allRatio_x,nanmedian(allRatio,1),...
        '--','Color',ones(1,3)*0.6,'markerfacecolor',ones(1,3)*0.7);
    ax = gca;
    ax.YColor = ones(1,3)*0.4;
    formatSubP2(hrat,hfut,dur)
    hold off
    yyaxis left
    formatSubP1(hpre,hfut,durVec([1]))
    % grid minor
    
    figPath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\ARF';
    fileName = sprintf('%s%sARF_NERC_JJA_%s_%02dh',figPath,filesep,regionName,dur);
    savePlot(fileName,'XYWH',[150,0,280,250],'needreply','N');
    
end
%% AUXILLARY FUNCTION
function colorMat = getColorMat(durVec,type)
switch(type)
    case 'pas'
        colorMat = [66 146 199]/255;
    case 'fut'
        colorMat = [240  60  43]/255;
end
end

function [A,B,arf] = getData(dataSourcePath,regionName,ENSEMBLENO,dur)
% unit: dur  [h]
region = getfield(REGIONS_info(),regionName);
[A,B] = deal([]);
for ensNo = ENSEMBLENO
    ensNo = ensNo{1};
    try
        fileName = sprintf('%s%sARFs_%s_ensNo%s_%s.mat',dataSourcePath,filesep,...
            region.Name,ensNo,'1980-2000');
        load(fileName,'T');try T.allT_Max = [];catch;end
        T(T.dur ~= dur,:) = [];
        T(~any(~isnan(T.arf'))',:) = [];
        A = [A;T];
        fileName = sprintf('%s%sARFs_%s_ensNo%s_%s.mat',dataSourcePath,filesep,...
            region.Name,ensNo,'2060-2080');
        load(fileName,'T');try T.allT_Max = [];catch;end
        T(T.dur ~= dur,:) = [];
        T(~any(~isnan(T.arf'))',:) = [];
        B = [B;T];
    catch me
        break
    end
end

arf = struct;
arf.area = [10,30,80,150,250,380,530,700,900,1250,2000,4000,8000]./(2.2).^2;
% arf.area = [10,30,80,150,250,380,530,700,900]./(2.2).^2;
end
function formatSubP2(hpre,hfut,durVec)
ylim([0.7,1.3])
grid minor
box on;
ax = gca;
ax.LineWidth = 1;
% ax.XScale = 'log';
ax.XTick = [100,1000,2000,3000,4000,5000];
ax.XLim = [10,5000];
plot(ax.XLim,[1,1],'k:','color',[0.5,0.5,0.5,0.5]);
% legend(flip([hpre]),flip(strcat(string(durVec'),'h')),'Location','SouthWest')
legend boxoff
xlabel('Area (km^2)')
ylabel('Future/Past');
end

function formatSubP1(hpre,hfut,durVec)
ylim([0.3,1])
grid minor
box on;
ax = gca;
legend(([hpre,hfut]),({strcat(string(durVec'),'h-past'),...
    strcat(string(durVec'),'h-future')}),...
    'Location','NorthEast')
legend boxoff
xlabel('Area (km^2)')
ylabel('Areal Reduction Factor');
% ax.XScale = 'log';
ax.XTick = [100,1000,2000,3000,4000,5000];
ax.XLim = [10,5000];
ax.LineWidth = 1;
end

