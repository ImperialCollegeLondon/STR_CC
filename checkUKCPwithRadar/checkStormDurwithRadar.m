
% colm = pink(8);
% colm(1:2,:) = [];
% colOBS = copper(8);
% colOBS(1:2,:) = [];
% lineType = 'x';%'-';
% lineTypeObs = '^';% 'o-';
% linewidth = 1;
% linewidthObs = 2;
% Name = {'London','SWestuk','Westuk','Scotland'};%{'London','SWestuk','Westuk','Scotland'};
dataSP = 'H:\DATA_CLIMATE\UKCP18\';

REGIONS = REGIONS_info();


Name = {'London','SWestuk','Westuk','Scotland'};

figure;
ha = tight_subplot(4,length(Name),[.08 .02],[.05 .05],[.10 .05]);
setFigureProperty('Paper')
set(0,'defaultAxesFontSize',14,'defaultAxesFontName','Arial',...
    'defaultAxesTitleFontWeight','Bold',...
    'defaultTextFontSize',10)

global seasons yi
seasons = 1:4;%1:4;%1:4
for thisstation = 1:length(Name)
    
    CPM = load(sprintf('%sDUR_Area_%s_1980-2000_4Season.mat',dataSP,Name{thisstation}),'dataplot');
    RAD = load(sprintf('%sObs_DUR_Area_%s_4Season.mat',dataSP,Name{thisstation}),'dataplot');
    
    fpIMF = 5;
    
    % PWAR
    yi = 0.01:0.01:1;% 0:0.001:1;
    axes(ha(thisstation))
    ax = plotOne(CPM,RAD,'pWAR',fpIMF);
    ax.YLim=[1e-1,6];%[0,1];%
    formatIt('Peak Wet Area Ratio',Name{thisstation})
    % titleIt(Name{thisstation});
    
    % median WAR
    yi = 0.01:0.01:1;% 0:0.001:1;
    axes(ha(thisstation+length(Name)));
    ax = plotOne(CPM,RAD,'WAR',fpIMF);
    ax.YLim=[1e-3,10];%[0,1];%
    formatIt('Q50 Wet Area Ratio',Name{thisstation})
    % title(Name{thisstation})
    
    % Dur
    yi = 1:1:48;%2:1:40;%
    axes(ha(thisstation+length(Name)*2));
    [ax,cpmNum,radNum] = plotOne(CPM,RAD,'Dur',fpIMF);
    ax.YLim = [3e-4,0.1];%[0,1];%
    ax.XLim = [1.5,48];
    formatIt('Duration',Name{thisstation})
    ax.XTick= [2,12,24,36,48];
    % title(Name{thisstation})
    
    % pIMF
    yi = 5:0.5:100;%5:1:100;%
    axes(ha(thisstation+length(Name)*3))
    ax = plotOne(CPM,RAD,'pIMF',fpIMF);
    ax.YLim=[1e-5,0.12];%[0,1];%
    formatIt('Peak Intensity',Name{thisstation})
    ax.XTick= [0,25,50,75,100];
    % title(Name{thisstation})
    
    %{
        % TotalRain
        for ikkk = 1
        yi = log2(0.01:0.1:10);%0.1:0.01:6;% log2(0.1:0.1:10);
        subplot(5,length(Name),thisstation+length(Name)*4)
        f = [];
        for ensNo = 1:12
            A = CPM.dataplot{1, thisseason}{ensNo}.IMF.*CPM.dataplot{1, thisseason}{ensNo}.Dur;
            [N,EDGES] = histcounts(A,[yi,inf],'Normalization','pdf');
            f(ensNo,:) = N;%
            % ksdensity(CPM.dataplot{1, thisseason}{ensNo}.IMF.*CPM.dataplot{1, thisseason}{ensNo}.Dur,...
            %     yi,'function','pdf','BoundaryCorrection','reflection');hold on;
        end
        plot(yi,log2(nanmedian(f,1)),'r-');
        ylim([log2(0.008),log2(2)])
        ax = gca;
        ax.YLim=log2([0.008,1.02]);
        text(ax.XLim(2),ax.YLim(2),'Total R[mm] ','HorizontalAlignment',...
            'Right','VerticalAlignment','Top');
        hsimran = fill([yi,flip(yi)],log2([nanmin(f,[],1),flip(nanmax(f,[],1))]),'r',...
            'LineStyle','none');hold on;
        
%         f=[];
%         for ensNo = 1:3
%             A = CPM_fut.dataplot{1, thisseason}{ensNo}.IMF.*CPM_fut.dataplot{1, thisseason}{ensNo}.Dur;
%             [N,EDGES] = histcounts(A,[yi,inf],'Normalization','pdf');
%             f(ensNo,:) = N;%
%             % ksdensity(CPM_fut.dataplot{1, thisseason}{ensNo}.IMF.*CPM_fut.dataplot{1, thisseason}{ensNo}.Dur,...
%             %     yi,'function','pdf','BoundaryCorrection','reflection');hold on;
%         end
%         plot(yi,log2(nanmedian(f,1)),'b-');
%         hsimfutran = fill([yi,flip(yi)],log2([nanmin(f,[],1),flip(nanmax(f,[],1))]),'b',...
%             'LineStyle','none');hold on;
        
        alpha(0.3)
        fobs = histcounts(RAD.dataplot{1, thisseason}.IMF.*RAD.dataplot{1, thisseason}.Dur,...
            [yi,inf],'Normalization','pdf');
        % ksdensity(RAD.dataplot{1, thisseason}.IMF.*RAD.dataplot{1, thisseason}.Dur,yi,'function','pdf',...
        %     'BoundaryCorrection','reflection');
        hobs = plot(yi,log2(fobs),'k-');
        box on;
        ylabel('PDF');
        set(ax,'linewidth',2)
        ax.YTick = log2(0.008):1:log2(2);
        ax.YTickLabel = round(1000*2.^(ax.YTick))/1000;
        xlim([0,6])
        end
    %}
    
    fprintf('Identified Events:\t %s CPM:%d/yr \t RAD:%d/yr \n',...
        Name{thisstation},round(nanmedian(cpmNum/20)),round(radNum/12));
    
    
end
% legend([hsimran,hsimfutran,hobs],{'CPM(1980-2000)','CPM(2060-2080)','Radar(2007-2018)'})

% legend([hsimran,hobs],{'CPM(1980-2000)','Radar(2007-2018)'})

rowno = 4;
colno = numel(Name);
set(ha(((1:rowno)-1)'*colno+(2:colno)),'YTickLabel',[],'YLabel',[])
% set(ha(1:(colno*(rowno-1))),'XTickLabel',[],'XLabel',[])

filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\';
filename = [filePath,filesep,'SeveralPDFs_event_',getSeasonName(seasons(1)),'to',getSeasonName(seasons(end))];
savePlot(filename,'units','centimeters','XYWH',[5,0,20,22],'needreply','Y');



function formatIt(xlab,name0)
ylabel('PDF');
set(gca,'linewidth',2)
% set(gca,'YScale','log')
% text(ax.XLim(2),ax.YLim(2),xlab,'HorizontalAlignment',...
%         'Right','VerticalAlignment','Top');
xlabel(xlab)
ax = gca;
text(ax.XLim(2)*0.9,ax.YLim(2),name0,'HorizontalAlignment',...
    'Right','VerticalAlignment','Top');
end

function titleIt(name0)
if strcmp(name0,'London')
    name = 'SouthUK';
end
title(sprintf('%s',name0))% (110x110Km^2)
end

function ax = plotOne_emp(CPM,RAD,var,fpIMF)
global seasons yi

sm = @(x)reshape(smooth(x,5),[1,numel(x)]);
hold on;
f = [];
for ensNo = 1:12
    A = [];
    for seai = reshape(seasons,1,[])
        A0 = getfield(CPM.dataplot{1,seai}{ensNo},var);
        rmind = CPM.dataplot{1,seai}{ensNo}.pIMF<fpIMF;
        A0(rmind) = [];
        A = cat(1,A,A0);
    end
    [N,EDGES] = histcounts(A,[yi,inf],'Normalization','pdf');
    f(ensNo,:) = N;% ksdensity(A,yi,'function','pdf');hold on;
end
ax = gca;
plot(yi,sm(nanmedian(f,1)),'r-');% ylim([0,0.25])
hsimran = fill([yi,flip(yi)],[sm(nanmin(f,[],1)),sm(flip(nanmax(f,[],1)))],'r',...
    'LineStyle','none');

alpha(0.3)
A = [];
for seai = reshape(seasons,1,[])
    A0 = getfield(RAD.dataplot{1,seai},var);
    rmind = RAD.dataplot{1,seai}.pIMF<fpIMF;
    A0(rmind) = [];
    A = cat(1,A,A0);
end
[N,EDGES] = histcounts(A,[yi,inf],'Normalization','pdf');
fobs = N;%ksdensity(getfield(RAD.dataplot{1,thisseason},var),yi,'function','pdf');
plot(yi,sm(fobs),'k-');% ,'linewidth',0.5);
box on;
hold off;

set(gca,'XScale','log')
end




function [ax,cpmNum,radNum] = plotOne(CPM,RAD,var,fpIMF)
global seasons yi
bw = 10;%%(yi(2)-yi(1))*10;
if yi(2)-yi(1) == 1
    bw = 1;
end
sm = @(x)x;%reshape(smooth(x,bw),size(x));
kernelName = 'normal';%'box';%

% hold on;

f = [];
cpmNum = [];
for ensNo = 1:12
    A = [];
    for seai = reshape(seasons,1,[])
        A0 = getfield(CPM.dataplot{1,seai}{ensNo},var);
        rmind = CPM.dataplot{1,seai}{ensNo}.pIMF<fpIMF;
        A0(rmind) = [];
        A = cat(1,A,A0);
    end
    cpmNum(ensNo) = numel(A);
    if yi(2)-yi(1) == 1
    f(ensNo,:) = ksdensity(A,yi,'function','pdf','kernel',kernelName,...
        'support', 'positive', 'BoundaryCorrection',...
            'Reflection','Bandwidth',2);
    else
        [f(ensNo,:),~,bw0]= ksdensity(A,yi,'function','pdf','kernel',kernelName,...
        'support', 'positive', 'BoundaryCorrection',...
            'Reflection');
        [f(ensNo,:),~,~]= ksdensity(A,yi,'function','pdf','kernel',kernelName,...
        'support', 'positive', 'BoundaryCorrection',...
            'Reflection','Bandwidth',bw0/2);
    end
end


% if yi(2)-yi(1) == 1
% 1
% end

ax = gca;
plot(yi,sm(nanmedian(f,1)),'r-');% ylim([0,0.25])
hold on;
hsimran = fill([yi,flip(yi)],[sm(nanmin(f,[],1)),sm(flip(nanmax(f,[],1)))],'r',...
    'LineStyle','none');
hold on;
alpha(0.3)
A = [];
radNum = [];
for seai = reshape(seasons,1,[])
    A0 = getfield(RAD.dataplot{1,seai},var);
    rmind = RAD.dataplot{1,seai}.pIMF<fpIMF;
    A0(rmind) = [];
    A = cat(1,A,A0);
end
radNum = numel(A);
if yi(2)-yi(1) == 1
    fobs = ksdensity(A,yi,'function','pdf',...
        'kernel',kernelName,'Bandwidth',2);
else
    [fobs,~,bw0] = ksdensity(A,yi,'function','pdf','kernel',kernelName);
    [fobs,~,~] = ksdensity(A,yi,'function','pdf','kernel',kernelName,...
        'support', 'positive', 'BoundaryCorrection',...
            'Reflection','Bandwidth',bw0/2);
end
plot(yi,sm(fobs),'k-');
box on;

hold off;



if strcmp(var,'dur')
    1
end

end


