% ----------------------------------------------------------------------- %
% THIS SCRIPT IS TO VISUALIZE SEVERAL STATISTICS FOR CONVECTIVE STROM DURING
% SUMMERTIME
%
% # data #
% ukcp18 cpm2.2 over the UK, all ensemble members
%
% # Conditional Density Plot #
%
% # storm extraction
% Current Criteria to extract convective storm is :
% all snapshots having Pmax>10mm/h will be taken as Convective Hour and used
% to plot figures below.
%
% # for each hour
% The case several hours from one convective storms are counted and
% included might happen.
% That's due to the statement by (Prein 2017) in # Nature # (if my
% understanding is correct). Paper Fig.2: he said all 'MCSs hours' are used.
% 
% # Study area
% When computing statistics for each sub-regions in the UK: only land area 
% is taken into account.
%
% @ Yuting Chen
% Imperial College London
% yuting.chen17@imperial.ac.uk
% update: 2020.08.06
% ----------------------------------------------------------------------- %

clear;clc
warning off

global regionName ensNo in

ENSEMBLENO = getEnsNos();

for regionName = {'CPM_NW','CPM_S','CPM_NE'}%{'WAL','EUK','SCO'}
    regionName = regionName{1};
    setFigureProperty('Paper');
    
    [E1,N1] = getEN(getfield(REGIONS_info(),regionName));
    UKMap = getUKMap();
    in = inpolygon(E1,N1,UKMap.borderE/1000,UKMap.borderN/1000);

    % Several Config
    for y_name = {'rami'}%{'rpmax'}%
        y_name = y_name{1};
        for FP = {'2060-2080'}%{'2020-2040','2060-2080'}%{'2007-2018'}%
            FP = FP{1};
            for x_name = {'rpmax','rsize','rspeed'}%{'rsize'}%
                x_name = x_name{1};
                
                [STATS0,STATS] = deal([]);
                for enInd = 1:length(ENSEMBLENO)
                    ensNo = ENSEMBLENO{enInd};
                    STATS0 = [STATS0;get4Plot('1980-2000',ensNo)];
                    STATS = [STATS;get4Plot(FP,ensNo)];
                end
                fprintf('%s\tPast%03.0f\tFuture%03.0f\t%2.1f\n',...
                    regionName,size(STATS,1)/12/20,...
                    size(STATS0,1)/12/20,...
                    100*(size(STATS0,1)-size(STATS,1))./(size(STATS,1)));
                
                [conpd1,samNum1,conpd2,samNum2] = deal([]);
                for enInd = 1:length(ENSEMBLENO)
                    ensNo = ENSEMBLENO{enInd};
                    thisSTATS = STATS0(strcmp(STATS0.ensNo,ensNo),:);
                    conX = getfield(thisSTATS,y_name);
                    Y = getfield(thisSTATS,x_name);
                    [X,Y,conpd1(enInd,:,:),samNum1(enInd,:,:)] = compute_CondiP(Y,conX,x_name,y_name);
                    
                    thisSTATS = STATS(strcmp(STATS.ensNo,ensNo),:);
                    conX = getfield(thisSTATS,y_name);
                    Y = getfield(thisSTATS,x_name);
                    [X,Y,conpd2(enInd,:,:),samNum2(enInd,:,:)] = compute_CondiP(Y,conX,x_name,y_name);
                end
                
                sampelNoThre = 12*5;
                
                [isSigChange] = testSig(conpd1,conpd2,samNum1,samNum2);
                reldif_conpd = 100*(conpd2-conpd1);
                reldif_conpd = squeeze(nanmedian(reldif_conpd,1))./squeeze(nanmedian(conpd1,1));
                samNum1 = squeeze(nansum(samNum1,1));
                samNum2 = squeeze(nansum(samNum2,1));
                
                [h] = plot_CondiChange(getfield(STATS0,x_name),getfield(STATS0,y_name),...
                    getfield(STATS,x_name),getfield(STATS,y_name),...
                    X,Y,reldif_conpd,samNum1,samNum2,sampelNoThre,isSigChange,...
                    x_name,y_name,FP);
                
                savePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\StormTracking';
                fileName = sprintf('UKCP_plot_CStorm_JJA_%s_con%s_%s_dif%s_%s-%s',x_name,y_name,regionName,...
                    FP,ENSEMBLENO{1},ENSEMBLENO{end});
                savePlot([savePath,filesep,fileName],'XYWH',[50,0,350,500],'needreply','N','onlyPng',true);
                
                pause(0.1)
                close all
            end
        end
    end
end


%%
function [isSigChange] = testSig(A,B,samNum1,samNum2)
% 1: different
% 0: not different.

% A(isnan(A)) = 0;
% B(isnan(B)) = 0;
c=A-B;
isSigChange = squeeze(ttest(A,B,'alpha',0.05));

end

%%
function [h] = plot_CondiChange(hist_x,hist_y,future_x,future_y,...
    PX0,PY0,reldif_conpd,samNum1,samNum2,sampelNoThre,isSigChange,...
    x_name,y_name,FP)

pascol = 'k';
futcol = [224	102	102]/255;

PMAP0 = reldif_conpd;% PMAP = conv2(PMAP,ones(5,5)/25,'same');

h = scatterhist(hist_x,hist_y,'Location','NorthEast',...
    'Direction','out','Color',pascol,...
    'LineWidth',[2,2,2]);
setXYAxis(h);
setBinwidth(h);

axes(h(1)); hold on; axis('square');
[i_s5,j_s5] = find(samNum2>sampelNoThre & isnan(PMAP0));
% plot(PX0(j_s5),PY0(i_s5),'color',pascol,'marker','.','linewidth',1,'markersize',3);

isSigChange = isSigChange*0+1;
[ij_ns] = find(isSigChange == 0 & ~isnan(PMAP0) & ...
    samNum2>sampelNoThre & samNum1>sampelNoThre);
PMAP0(ij_ns) = 0;

samNum1 = imresize(samNum1,'Scale',10,'Method','bilinear');
samNum2 = imresize(samNum2,'Scale',10,'Method','bilinear');
PMAP = imresize(PMAP0,'Scale',10,'Method','bilinear');
PX = linspace(PX0(1),PX0(end),size(PMAP,2));
PY = linspace(PY0(1),PY0(end),size(PMAP,1));
PMAP(samNum1<sampelNoThre | samNum2<sampelNoThre) = NaN;

[XX,YY] = meshgrid(PX,PY);
XX = XX(samNum2>sampelNoThre);
YY = YY(samNum2>sampelNoThre);
K = boundary(XX,YY,0.8);
plot(XX(K),YY(K),'--','color',pascol,'linewidth',1)
pcolor(PX,PY,PMAP);
formatContour();

c = formatColormap();
plotHist(h);
formatHist(h);
formatPosition(h);

    function setXYAxis(h)
        h(1).YAxisLocation = 'left';
        h(1).XAxisLocation = 'bottom';
        h(2).Children(1).EdgeColor = 'none';
        h(3).Children(1).EdgeColor = 'none';
    end

    function setBinwidth(h)
        if strcmp(x_name,'rspeed')
            h(2).Children(1).BinWidth = 5;
        elseif strcmp(x_name,'rsize')
            h(2).Children(1).BinWidth = 500;
        elseif strcmp(x_name,'rpmax')
            h(2).Children(1).BinWidth = 4;
        else
            % h(2).Children(1).BinWidth = h(2).Children(1).BinWidth*2;
        end
        h(3).Children(1).BinWidth = h(3).Children(1).BinWidth*5;% 2000;% h(3).Children(1).BinWidth*2;
        h(1).Children.Marker = 'none';
    end
    function formatContour()
        shading flat
        set(gca,'XLim',[PX(1),PX(end)])
        set(gca,'YLim',[PY(1),PY(end)])
        
        xlabel(getLabel(x_name))
        ylabel(getLabel(y_name))
        axis on
    end
    function plotHist(h)
        axes(h(2));hold on;% 'FaceColor',futcol,'DisplayStyle','stairs','EdgeColor',futcol
        histogram(future_x,'Normalization','pdf',...
            'BinWidth',h(2).Children(1).BinWidth,...
            'linewidth',1,'DisplayStyle','stairs','EdgeColor',futcol);
        alpha(0.4)
        
        axes(h(3));hold on;%'FaceColor',futcol,'DisplayStyle','stairs','EdgeColor',futcol
        histogram(future_y,'Normalization','pdf',...
            'BinWidth',h(3).Children(1).BinWidth,...
            'linewidth',1,'DisplayStyle','stairs','EdgeColor',futcol);
        alpha(0.4)
    end
    function c = formatColormap()
        cptcmap('diff_darkBlue_darkRed','mapping','scaled','ncol',21);
        c = colorbar('Location','SouthOutside');
        c.Ticks = [-400,-200,0,200,400];caxis([-420,420]);% #MIGHT NEED CHANGE#
        % c.Ticks = [-200,-050,0,050,200];caxis([-210,210])
        c.TickLabels = strcat(c.TickLabels,'%');
        c.TickLabels{3} = '0';
        set(gca,'TickDir','out','Linewidth',1,'Fontsize',14);
    end
    function formatHist(h)
        legend(h(2),[h(2).Children],{FP,'1980-2000'},'location','NorthOutside',...
            'Orientation','horizontal')
        axes(h(2));legend boxoff
        h(3).XLim = h(1).YLim; h(3).YLim = h(3).YLim*1.2;
        h(2).XLim = h(1).XLim; h(2).YLim = h(2).YLim*1.2;
        % axis(h(1),'auto');  % Sync axes
        set(gcf,'units','points','position',[50,0,350,500]);
        
        set(h(3),'YScale','log')
        set(h(2),'YScale','log')
    end
    function formatPosition(h)
        h(1).Position(2) = 0.25;
        if ~ strcmp(x_name,'rsize')
            h(1).Position(1) = h(1).Position(1)+0.05;
            h(3).Children(2).Visible = 'off';
            h(3).Children(1).Visible = 'off';
            h(1).Position(4) = h(1).Position(3)/500*350;
        else
            h(1).Position(1) = h(1).Position(1)+0.05;
            h(1).Position(4) = h(1).Position(3)/500*350;
        end
    end
    function label = getLabel(x_name)
        switch(x_name)
            case 'rdur'
                label = 'Duration(h)';
            case 'rspeed'
                label = 'Speed(km h^{-1})';
            case 'rvol'
                label = 'Ptotal(m^3 s^{-1})';
            case 'rpmax'
                label = 'Pmax(mm h^{-1})';
            case 'rsize'
                label = 'Size(km^2)';
            case 'rsizeall'
                label = 'Wet Area(km^2)';
            case 'rrmi'
                label = 'Areal Mean Rain Rate(mm)';
            case 'rami'
                label = 'Areal Mean Intensity(mm/h)';
            otherwise
                error('check Function getLabel(x_name)');
        end
        
    end
end

% AUXILLARY FUNCTION

function [X,Y,conpd,samNum] = compute_CondiP(Y,conX,varNX,varNY)

[gridx1,yWin] = getIt(varNX);
[gridx2,conWin] = getIt(varNY);

[X,Y,conpd,samNum] = smoothConDist([Y(:),conX(:)],...
    gridx1,gridx2,conWin,yWin);
conpd(isnan(conpd))=0;

    function [gridx,winx]=getIt(varNX)
        switch(varNX)
            case 'rvol'
                gridx = 0:800:100000;
                winx = 8000;
            case 'rsize'
                gridx = 0:200:18000;%18000
                winx = 2000;
            case 'rsizeall'
                gridx = 0:300:100000;
                winx = 3000;
            case 'rspeed'
                gridx = 0:2:160;
                winx = 20;
            case 'rpmax'
                gridx = 10:2:130;
                winx = 20;
            case 'rdur'
                gridx = 1:2:24*7;
                winx = 20;
            case 'rrmi'
                gridx = 1:0.1:10;
                winx = 0.5;
            case 'rami'
                gridx = 0:0.02:2.5;
                winx = 0.2;
            otherwise
                error(sprintf('Please add info for:%s',varname));
        end
    end

    function [X,Y,GMAP,samNum] = smoothConDist(XY,gridx1,gridx2,conWin,yWin)
        
        [samNum,hxy_y] = deal([]);
        
        % figure;
        x = XY(:,1);y = XY(:,2);
        
        tag = 1;
        for y_cur = gridx2
            y_cur_range = [y_cur - conWin/2,y_cur + conWin/2];
            x_cur = x(y < y_cur_range(2) & y > y_cur_range(1));
            % smooth window
            x_cur = histcounts(x_cur,[-inf,(gridx1(1:end-1)+gridx1(2:end))/2,inf]);
            conVec = round(yWin/(gridx1(2)-gridx1(1)));
            if mod(conVec,2)==0
                conVec = conVec+1;
            end
            x_cur = conv(ones(1,conVec),x_cur);
            x_cur = x_cur(((conVec-1)/2+1):end-((conVec-1)/2));
            
            samNum(tag,:) = x_cur(:);
            hxy_y(tag,:) = x_cur(:)./nansum(x_cur);
            tag = tag+1;
        end
        
        PMAP = hxy_y;
        PMAP(PMAP==0) = NaN;
        X = gridx1;
        Y = gridx2;
        GMAP = PMAP;
    end
end

function [A] = get4Plot(Period,ensNo)
global regionName in

STATS = [];
MON = 6:8;
Config = getConfig(upper(regionName),MON,Period,ensNo);
if strcmp(Period,'2007-2018')
    Dat = load([Config.saveIt.path,filesep,sprintf('CS_%s_%s_STATS_%02d-%02d_%s.mat',...
        regionName,Period,Config.Month(1),Config.Month(end),'RAD')],...
        'STATS','Config');
else
    Dat = load([Config.saveIt.path,filesep,sprintf('CS_%s_%s_STATS_%02d-%02d_%s.mat',...
        regionName,Period,Config.Month(1),Config.Month(end),ensNo)],...
        'STATS','Config');
end
STATS = [STATS;Dat.STATS];


% A = STATS(STATS.rpmax>=5 & ~isnan(STATS.rspeed),:);
% A = STATS(STATS.rpmax>=5 & ~isnan(STATS.rspeed),:);
A = STATS(STATS.rpmax>=10,:);
A.rami = A.rvol /getSumArea(regionName,in) /1000 *3600; % mm/h

% % STATS = STATS(STATS.rpmax>=5,:);
% % % STATS.tag = STATS.evi;
% % STATS.mon = [];
% % D = grpstats(STATS,{'evi'},{'mean','max','median'});
% % A = table();
% % A.rvol = D.max_rvol; % mean_rvol.*D.median_rdur/24;%
% % A.rsize = D.max_rsize;
% % A.rpmax = D.max_rpmax;
% % A.rspeed = D.median_rspeed;
% % A.rdur = D.median_rdur;

A.ensNo = repmat(string(ensNo),size(A.rpmax));

end
