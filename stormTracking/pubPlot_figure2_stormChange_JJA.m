% ----------------------------------------------------------------------- %
% THIS SCRIPT IS TO VISUALIZE SEVERAL STATISTICS FOR CONVECTIVE STROM DURING
% SUMMERTIME
%
% # some necessary description #
%
% # Conditional Density Plot #
%
% # storm extraction
% Current Criteria to extract convective storm is :
% all snapshots having Pmax>5mm/h will be taken as Convective Hour and used
% to plot figures below.
% The reason of that is this threshold were used in several cell tracking
% algorithm for example in storm tracking using radar (35 dBZ).
%
% # for each hour
% The case several hours from one convective storms are counted and
% included might happen.
% That's due to the statement by (Prein 2017) in # Nature # (if my
% understanding is correct). Paper Fig.2: he said all 'MCSs hours' are used.
%
%
% @ Yuting Chen
% Imperial College London
% yuting.chen17@imperial.ac.uk
% ----------------------------------------------------------------------- %
clear;clc
warning off
global regionName ensNo
ENSEMBLENO = getEnsNos();
for regionName = {'WAL','EUK','SCO'}
    regionName = regionName{1};
    setFigureProperty('Paper');
    % Several Config
    for y_name = {'rvol'};%{'rpmax'}% {'rsize'}%{'rdur'}%  {'rdur','rvol'};
        y_name = y_name{1};
        for FP = {'2060-2080'}%{'2020-2040','2060-2080'}%{'2007-2018'}%
            FP = FP{1};
            % figure;
            % hp_rpmax = uipanel('position',[0 0 0.25 1]);
            % hp_rsize = uipanel('position',[0.25 0 0.25 1]);
            % hp_rspeed = uipanel('position',[0.5 0 0.25 1]);
            % hp_rvol = uipanel('position',[0.75 0 0.25 1]);
            for x_name = {'rdur','rvol','rpmax','rsize','rspeed'}% {'rpmax','rsize','rspeed'}%,'rdur'}%,'rdur'}
                x_name = x_name{1};
                % eval(['hp = hp_',x_name,';']);
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
                    [X,Y,conpd1(enInd,:,:),samNum1(enInd,:,:)] = get_CondiP(Y,conX,x_name,y_name);
                    
                    thisSTATS = STATS(strcmp(STATS.ensNo,ensNo),:);
                    conX = getfield(thisSTATS,y_name);
                    Y = getfield(thisSTATS,x_name);
                    [X,Y,conpd2(enInd,:,:),samNum2(enInd,:,:)] = get_CondiP(Y,conX,x_name,y_name);
                end
                %%
                sampelNoThre = 12*5;
                
                [isSigChange] = testSig(conpd1,conpd2,samNum1,samNum2);
                reldif_conpd = 100*(conpd2-conpd1);%./conpd1;
                reldif_conpd = squeeze(nanmedian(reldif_conpd,1))./squeeze(nanmedian(conpd1,1));
                samNum1 = squeeze(nansum(samNum1,1));
                samNum2 = squeeze(nansum(samNum2,1));
                
                [h] = plot_CondiChange(getfield(STATS0,x_name),getfield(STATS0,y_name),...
                    getfield(STATS,x_name),getfield(STATS,y_name),...
                    X,Y,reldif_conpd,samNum1,samNum2,sampelNoThre,isSigChange,...
                    x_name,y_name,FP);
                set(gcf,'units','points','position',[50,0,350,500]);

                h(1).Position(2) = 0.3;
                if ~ strcmp(x_name,'rspeed')% hp.Position(1)<0.75
                    h(1).Position(1) = h(1).Position(1)+0.05;
                    h(3).Children(2).Visible = 'off';
                    h(3).Children(1).Visible = 'off';
                    h(1).Position(4) = h(1).Position(3)/500*350;
                else
                    h(1).Position(1) = h(1).Position(1)+0.05;
                    h(1).Position(4) = h(1).Position(3)/500*350;
                end
                
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

function [h] = plot_CondiChange(hist_x,hist_y,future_x,future_y,...
    PX0,PY0,reldif_conpd,samNum1,samNum2,sampelNoThre,isSigChange,...
    x_name,y_name,FP)

pascol = 'k';%[59  76 192]/255;
futcol = [191   7  41]/255;

PMAP0 = reldif_conpd;% PMAP = conv2(PMAP,ones(5,5)/25,'same');

h = scatterhist(hist_x,hist_y,'Location','NorthEast',...
    'Direction','out','Color',pascol,...
    'LineWidth',[2,2,2]);
h(1).YAxisLocation = 'left';
h(1).XAxisLocation = 'bottom';
if strcmp(x_name,'rspeed')
    h(2).Children(1).BinWidth = 5;
elseif strcmp(x_name,'rsize')
    h(2).Children(1).BinWidth = 1000;
elseif strcmp(x_name,'rpmax')
    h(2).Children(1).BinWidth = 5;
else
    h(2).Children(1).BinWidth = h(2).Children(1).BinWidth*2;
end
h(3).Children(1).BinWidth = 2000;% h(3).Children(1).BinWidth*2;
h(1).Children.Marker = 'none';
axes(h(1));hold on;

[i_s5,j_s5] = find(samNum2>sampelNoThre & isnan(PMAP0));
plot(PX0(j_s5),PY0(i_s5),'color',pascol,'marker','.','linewidth',1,'markersize',3);


% [i_s5,j_s5] = find(samNum2>sampelNoThre & samNum1>sampelNoThre & ...
%     isSigChange == 0 & abs(PMAP0)>10);
% plot(PX0(j_s5),PY0(i_s5),'k.','linewidth',1,'markersize',4);

isSigChange = isSigChange*0+1;%%%%%%%
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

shading flat
set(gca,'XLim',[PX(1),PX(end)])
set(gca,'YLim',[PY(1),PY(end)])
% cptcmap('cool-warm','mapping','scaled','ncol',21);
cptcmap('diff_darkBlue_darkRed','mapping','scaled','ncol',21);
caxis([-420,420]);% #MIGHT NEED CHANGE#
% caxis([-105,105])

xlabel(getLabel(x_name))
ylabel(getLabel(y_name))
axis on
c = colorbar('Location','SouthOutside');
c.Ticks = [-400,-200,0,200,400];
% c.Ticks = [-100,-50,-20,0,20,50,100];%%%%%%%%%%%%
c.TickLabels = strcat(c.TickLabels,'%');
c.TickLabels{3} = '0';
% alpha(0.7)
set(gca,'TickDir','out','Linewidth',1,'Fontsize',14);

axes(h(2));hold on;
histogram(future_x,'FaceColor',futcol,'Normalization','pdf',...
    'BinWidth',h(2).Children(1).BinWidth);
alpha(0.2)

axes(h(3));hold on;
histogram(future_y,'FaceColor',futcol,'Normalization','pdf',...
    'BinWidth',h(3).Children(1).BinWidth);
set(h(3),'YScale','log')
set(h(2),'YScale','log')
alpha(0.2)

legend(h(2),[h(2).Children],{FP,'1980-2000'},'location','NorthOutside',...
    'Orientation','horizontal')
axes(h(2));legend boxoff
h(3).XLim = h(1).YLim; h(3).YLim = h(3).YLim*1.2;
h(2).XLim = h(1).XLim; h(2).YLim = h(2).YLim*1.2;
% axis(h(1),'auto');  % Sync axes
end

% AUXILLARY FUNCTION
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
    otherwise
        error('check Function getLabel(x_name)');
end

end

function [A] = get4Plot(Period,ensNo)
global regionName

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
A = STATS(STATS.rpmax>=5 & ~isnan(STATS.rspeed),:);


% STATS = STATS(STATS.rpmax>=5,:);
% % STATS.tag = STATS.evi;
% STATS.mon = [];
% D = grpstats(STATS,{'evi'},{'mean','max','median'});
% A = table();
% A.rvol = D.max_rvol; % mean_rvol.*D.median_rdur/24;%
% A.rsize = D.max_rsize;
% A.rpmax = D.max_rpmax;
% A.rspeed = D.median_rspeed;
% A.rdur = D.median_rdur;

A.ensNo = repmat(string(ensNo),size(A.rpmax));

end

function [X,Y,conpd,samNum] = get_CondiP(Y,conX,varNX,varNY)

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
                gridx = 0:300:30000;
                winx = 3000;
            case 'rspeed'
                gridx = 0:2:160;
                winx = 20;
            case 'rpmax'
                gridx = 0:2:130;
                winx = 20;
            case 'rdur'
                gridx = 1:2:24*7;
                winx = 20;
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
