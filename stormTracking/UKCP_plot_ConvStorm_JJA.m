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

global regionName ensNo
ENSEMBLENO = getEnsNos();
for regionName = {'WAL','EUK','SCO'}
    regionName = regionName{1};
    setFigureProperty('Paper');
    % Several Config
    for y_name = {'rdur'}% {'rdur','rvol'};
        y_name = y_name{1};
        for FP = {'2060-2080'}%{'2020-2040','2060-2080'}%{'2007-2018'}%
            FP = FP{1};
            % figure;
            % hp_rpmax = uipanel('position',[0 0 0.25 1]);
            % hp_rsize = uipanel('position',[0.25 0 0.25 1]);
            % hp_rspeed = uipanel('position',[0.5 0 0.25 1]);
            % hp_rvol = uipanel('position',[0.75 0 0.25 1]);
            for x_name = {'rpmax','rsize','rspeed','rvol'}%,'rdur'}
                x_name = x_name{1};
                % eval(['hp = hp_',x_name,';']);
                [STATS0,STATS] = deal([]);
                for ensNo = ENSEMBLENO
                    ensNo = ensNo{1};
                    STATS0 = [STATS0;get4Plot('1980-2000',ensNo)];
                    STATS = [STATS;get4Plot(FP,ensNo)];
                end
                PL = 0;
                % figure;
                conX = getfield(STATS0,y_name);
                Y = getfield(STATS0,x_name);
                [h,X,Y,conpd1,samNum1] = plot_CondiP(Y,conX,x_name,y_name,PL);
                % title(sprintf('%s-JJA','1980-2000'));
                % figure;
                conX = getfield(STATS,y_name);
                Y = getfield(STATS,x_name);
                [h,X,Y,conpd2,samNum2] = plot_CondiP(Y,conX,x_name,y_name,PL);
                % title(sprintf('%s-JJA',FP));
                
                
                sampelNoThre = 5;
                [h] = plot_CondiChange(getfield(STATS0,x_name),getfield(STATS0,y_name),...
                    getfield(STATS,x_name),getfield(STATS,y_name),...
                    X,Y,conpd1,samNum1,conpd2,samNum2,sampelNoThre,...
                    x_name,y_name,FP);
                set(gcf,'units','points','position',[50,0,350,500]);

                h(1).Position(2) = 0.3;
                if ~ strcmp(x_name,'rvol')% hp.Position(1)<0.75
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
                savePlot([savePath,filesep,fileName],'XYWH',[50,0,350,500],'needreply','N','onlyPng',false);
                pause(0.1)
                close all
            end
        end
    end
end


%%
function [h] = plot_CondiChange(hist_x,hist_y,future_x,future_y,...
    PX,PY,conpd1,samNum1,conpd2,samNum2,sampelNoThre,x_name,y_name,FP)

PMAP = 100*(conpd2-conpd1)./conpd1;% PMAP = conv2(PMAP,ones(5,5)/25,'same');
sampelNoThre = 5;
PMAP(samNum1<sampelNoThre | samNum2<sampelNoThre) = NaN;

period = [hist_x*0;future_x*1];
x = [hist_x;future_x];
y = [hist_y;future_y];

% this one is SUPER SLOW
% scatterhist(x,y,'Group',period,'Location','SouthEast',...
%     'Direction','out','Color','kb','LineStyle',{'-','-.'},...
%     'LineWidth',[2,2,2],'Marker','+od','MarkerSize',[4,5]);
% this one is SUPER SLOW

h = scatterhist(hist_x,hist_y,'Location','NorthEast',...
    'Direction','out','Color','kb',...
    'LineWidth',[2,2,2]);
h(1).YAxisLocation = 'left';
h(1).XAxisLocation = 'bottom';
if strcmp(x_name,'rspeed')
    h(2).Children(1).BinWidth = h(2).Children(1).BinWidth*4;
else
    h(2).Children(1).BinWidth = h(2).Children(1).BinWidth*2;
end
h(3).Children(1).BinWidth = 2;% h(3).Children(1).BinWidth*2;
h(1).Children.Marker = 'none';
axes(h(1));hold on;

[i_s5,j_s5] = find(samNum2>sampelNoThre & isnan(PMAP));
plot(PX(j_s5),PY(i_s5),'k.','linewidth',1,'markersize',2);hold on


samNum2 = imresize(samNum2,'Scale',[10,10],'Method','bilinear');
PMAP = imresize(PMAP,'Scale',[10,10],'Method','bilinear');
PX = linspace(PX(1),PX(end),size(PMAP,2));
PY = linspace(PY(1),PY(end),size(PMAP,1));

[XX,YY] = meshgrid(PX,PY);
XX = XX(samNum2>sampelNoThre);
YY = YY(samNum2>sampelNoThre);

K = boundary(XX,YY,0.5);
plot(XX(K),YY(K),'k--','linewidth',1)
pcolor(PX,PY,PMAP);
shading flat
set(gca,'XLim',[PX(1),PX(end)])
set(gca,'YLim',[PY(1),PY(end)])
cptcmap('diff_darkBlue_darkRed','mapping','scaled','ncol',11);
caxis([-220,220]);% #MIGHT NEED CHANGE#

xlabel(getLabel(x_name))
ylabel(getLabel(y_name))
axis on
c = colorbar('Location','SouthOutside');
c.Ticks = [-200,-100,0,100,200];
c.TickLabels = strcat(c.TickLabels,'%');
set(gca,'TickDir','out','Linewidth',1,'Fontsize',14);

axes(h(2));hold on;
histogram(future_x,'FaceColor','r','Normalization','pdf',...
    'BinWidth',h(2).Children(1).BinWidth);
alpha(0.2)

axes(h(3));hold on;
histogram(future_y,'FaceColor','r','Normalization','pdf',...
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
        label = 'Ptotal(mm^3 s^{-1})';
    case 'rpmax'
        label = 'Pmax(mm h^{-1})';
    case 'rsize'
        label = 'Size(km^2)';
    otherwise
        error('check Function getLabel(x_name)');
end

end

function A = get4Plot(Period,ensNo)
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

A = STATS(STATS.rpmax>=5 & ~isnan(STATS.rspeed),:);

A = table();
STATS.tag = STATS.evi;
STATS.mon = [];
D = grpstats(STATS,{'tag'},{'mean','max','median'});
A.rvol = D.max_rvol; % mean_rvol.*D.median_rdur/24;%
A.rsize = D.max_rsize;
A.rpmax = D.max_rpmax;
A.rspeed = D.median_rspeed;
A.rdur = D.median_rdur;

A.ensNo = repmat(string(ensNo),size(A.rpmax));

end

function [h,X,Y,conpd,samNum] = plot_CondiP(Y,conX,varNX,varNY,PL)

[gridx1,yWin] = getIt(varNX);
[gridx2,conWin] = getIt(varNY);

if PL
    h = scatterhist(Y,conX);
    h(1).Children.Marker = 'none';
    axes(h(1))
else
    h = [];
end
[X,Y,conpd,samNum] = plotJointDist([Y(:),conX(:)],...
    gridx1,gridx2,conWin,yWin,PL);
conpd(isnan(conpd))=0;
xlabel(varNX)
ylabel(varNY)
set(gca,'TickDir','out')

    function [gridx,winx]=getIt(varNX)
        switch(varNX)
            case 'rvol'
                gridx = 0:600:80000;
                winx = 3000;
            case 'rsize'
                gridx = 0:800:160000;
                winx = 5000;
            case 'rspeed'
                gridx = 0:2:50;
                winx = 10;
            case 'rpmax'
                gridx = 0:2:150;
                winx = 10;
            case 'rdur'
                gridx = 1:2:24*3;
                winx = 10;
            otherwise
                error(sprintf('Please add info for:%s',varname));
        end
    end
end

function [X,Y,GMAP,samNum] = plotJointDist(XY,gridx1,gridx2,conWin,yWin,pl)


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

% [x1,x2] = meshgrid(gridx1, gridx2);
% x1 = x1(:);
% x2 = x2(:);
% xi = [x1 x2];
% hxy = ksdensity([x,y],xi,'BoundaryCorrection','reflection','PlotFcn','surf','Kernel','epanechnikov');
% hxy = reshape(hxy,numel(gridx2),numel(gridx1));
%
% samNum = histcounts2(x,y,'XBinEdges',[-inf,(gridx1(1:end-1)+gridx1(2:end))/2,inf],...
%     'YBinEdges',[-inf,(gridx2(1:end-1)+gridx2(2:end))/2,inf]);
% hxy_y = hxy./reshape(nansum(hxy,2),[size(hxy,1),1]);
PMAP = hxy_y;

if pl==1
    % surf(gridx1,gridx2,PMAP);
    hold on
    pcolor(gridx1,gridx2,PMAP);shading interp;
    shading flat
    set(gca,'XLim',[min(x) max(x)])
    set(gca,'YLim',[min(y) max(y)])
    cptcmap('blue_continuous','mapping','scaled','ncol',100,'flip',true);
    caxis([0,0.2]);% #MIGHT NEED CHANGE#
    xlabel('..')
    ylabel('Precipitation Volume')
    % ax = gca;
    % ax.XTickLabel = ax.XTick;
    % ax.YTickLabel = ax.YTick;
    axis on
    box on
end
PMAP(PMAP==0) = NaN;
X = gridx1;
Y = gridx2;
GMAP = PMAP;
box on;

end

% function plotMarginalHist(hist_varname,yi,STATS,STATS0)
% % hist_varname = 'rspeed';%'rpmax';%'rsize';%'rvol';%
% % yi = 0:2:30;%0:3000:120000;
% % f = ksdensity(getfield(STATS,hist_varname),yi,'function','pdf','Bandwidth',300);hold on;
% % plot(yi,f,'r-');%
% hold on;
% [N,EDGES] = histcounts(getfield(STATS,hist_varname),yi,'Normalization','pdf');hold on;
% bar((EDGES(1:end-1)+EDGES(2:end))/2,N,'r');alpha(0.2)
%
% % f = ksdensity(getfield(STATS0,hist_varname),yi,'function','pdf','Bandwidth',300);hold on;
% % plot(yi,f,'k-');
% [N,EDGES] = histcounts(getfield(STATS0,hist_varname),yi,'Normalization','pdf');hold on;
% bar((EDGES(1:end-1)+EDGES(2:end))/2,N,'k');alpha(0.2)
% hold off;
% legend('future','1980-2000')%
% xlabel(hist_varname)
% set(gca,'TickDir','out');
% end
