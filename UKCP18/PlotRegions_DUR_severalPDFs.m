setFigureProperty('Paper')
figure;
Name = {'London','Westuk','Scotland'};%{'London','SWestuk','Westuk','Scotland'};
set(0,'defaultAxesFontSize',14,'defaultAxesFontName','Arial',...
    'defaultAxesTitleFontWeight','Bold',...
    'defaultTextFontSize',10)
global thisseason yi
for thisseason = 2%1:4
    for thisstation = 1:length(Name)
        
        CPM = load(sprintf('%sDUR_Area_%s_1980-2000_4Season.mat',dataSP,Name{thisstation}),'dataplot');
        RAD = load(sprintf('%sObs_DUR_Area_%s_4Season.mat',dataSP,Name{thisstation}),'dataplot');
        CPM_fut = load(sprintf('%sDUR_Area_%s_2060-2080_4Season.mat',dataSP,Name{thisstation}),'dataplot');
        
        fpIMF = 5;
        
        % PWAR
        yi = 0.01:0.01:1;% 0:0.001:1;
        subplot(4,length(Name),thisstation)
        ax = plotOne(CPM,CPM_fut,RAD,'pWAR',fpIMF);
        ax.YLim=[1e-1,6];%[0,1];%
        text(ax.XLim(2),ax.YLim(2),'Peak Wet Area Ratio','HorizontalAlignment',...
            'Right','VerticalAlignment','Top');
        ylabel('PDF');
        set(gca,'linewidth',2)
        set(gca,'YScale','log')
        if strcmp(Name{thisstation},'London')
            title(sprintf('%s (110x110Km^2)','SouthUK'))
        else
            title(sprintf('%s (110x110Km^2)',Name{thisstation}))
        end
        % median WAR
        yi = 0.01:0.01:1;% 0:0.001:1;
        subplot(4,length(Name),thisstation+length(Name))
        ax = plotOne(CPM,CPM_fut,RAD,'WAR',fpIMF);
        ax.YLim=[1e-3,10];%[0,1];%
        text(ax.XLim(2),ax.YLim(2),'Q50 Wet Area Ratio','HorizontalAlignment',...
            'Right','VerticalAlignment','Top');
        ylabel('PDF');
        set(gca,'linewidth',2)
        set(gca,'YScale','log')
        % title(Name{thisstation})
        
        % Dur
        yi = 2:1:40;%2:1:40;% 
        subplot(4,length(Name),thisstation+length(Name)*2)
        ax = plotOne(CPM,CPM_fut,RAD,'Dur',fpIMF);
        ax.YLim = [3e-4,0.3];%[0,1];%
        ax.XLim = [2,40];
        text(ax.XLim(2),ax.YLim(2),'Duration','HorizontalAlignment',...
            'Right','VerticalAlignment','Top');
        ylabel('PDF');
        set(gca,'linewidth',2)
        set(gca,'YScale','log')
        % title(Name{thisstation})
        
        % pIMF
        yi = 5:0.5:100;%5:1:100;%
        subplot(4,length(Name),thisstation+length(Name)*3)
        ax = plotOne(CPM,CPM_fut,RAD,'pIMF',fpIMF);
        ax.YLim=[1e-5,0.15];%[0,1];%
        text(ax.XLim(2),ax.YLim(2),'Peak Intensity','HorizontalAlignment',...
            'Right','VerticalAlignment','Top');
        ylabel('PDF');
        set(gca,'linewidth',2)
        set(gca,'YScale','log')
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
    end
    % legend([hsimran,hsimfutran,hobs],{'CPM(1980-2000)','CPM(2060-2080)','Radar(2007-2018)'})
    
    % legend([hsimran,hobs],{'CPM(1980-2000)','Radar(2007-2018)'})
    
    
    filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\';
    filename = [filePath,filesep,'SeveralPDFs_event_',getSeasonName(thisseason)];
    savePlot(filename,'XYWH',[150,0,1200,800],'needreply','N');
    pause(2)
    % close all
    
end


function ax = plotOne_emp(CPM,CPM_fut,RAD,var,fpIMF)
global thisseason yi

sm = @(x)reshape(smooth(x,1),[1,numel(x)]);
hold on;
f = [];
for ensNo = 1:12
    A = getfield(CPM.dataplot{1,thisseason}{ensNo},var);
    rmind = CPM.dataplot{1,thisseason}{ensNo}.pIMF<fpIMF;
    A(rmind) = [];
    [N,EDGES] = histcounts(A,[yi,inf],'Normalization','pdf');
    f(ensNo,:) = N;% ksdensity(A,yi,'function','pdf');hold on;
end
ax = gca;
plot(yi,sm(nanmedian(f,1)),'r-');% ylim([0,0.25])
hsimran = fill([yi,flip(yi)],[sm(nanmin(f,[],1)),sm(flip(nanmax(f,[],1)))],'r',...
    'LineStyle','none');

f=[];
for ensNo = 1:3
    A = getfield(CPM_fut.dataplot{1,thisseason}{ensNo},var);
    rmind = CPM_fut.dataplot{1,thisseason}{ensNo}.pIMF<fpIMF;
    A(rmind) = [];
    [N,EDGES] = histcounts(A,[yi,inf],'Normalization','pdf');
    f(ensNo,:) = N;% ksdensity(A,yi,'function','pdf');hold on;
end
plot(yi,sm(nanmedian(f,1)),'b-');
hsimran = fill([yi,flip(yi)],[sm(nanmin(f,[],1)),sm(flip(nanmax(f,[],1)))],'b',...
    'LineStyle','none');


alpha(0.3)

A = getfield(RAD.dataplot{1,thisseason},var);
rmind = RAD.dataplot{1,thisseason}.pIMF<fpIMF;
A(rmind) = [];
[N,EDGES] = histcounts(A,[yi,inf],'Normalization','pdf');
fobs = N;%ksdensity(getfield(RAD.dataplot{1,thisseason},var),yi,'function','pdf');
plot(yi,sm(fobs),'k-');% ,'linewidth',0.5);
box on;
hold off;

set(gca,'XScale','log')
end




function ax = plotOne(CPM,CPM_fut,RAD,var,fpIMF)
global thisseason yi
bw = 10;%%(yi(2)-yi(1))*10;
if yi(2)-yi(1) == 1
    bw = 1;
end
sm = @(x)x%reshape(smooth(x,bw),size(x));
kernelName = 'normal';%'box';%

hold on;

f = [];
for ensNo = 1:12
    A = getfield(CPM.dataplot{1,thisseason}{ensNo},var);
    rmind = CPM.dataplot{1,thisseason}{ensNo}.pIMF<fpIMF;
    A(rmind) = [];
    f(ensNo,:) = ksdensity(A,yi,'function','pdf','kernel',kernelName);
end
ax = gca;
plot(yi,sm(nanmedian(f,1)),'r-');% ylim([0,0.25])

hsimran = fill([yi,flip(yi)],[sm(nanmin(f,[],1)),sm(flip(nanmax(f,[],1)))],'r',...
    'LineStyle','none');

f=[];
for ensNo = 1:3
    A = getfield(CPM_fut.dataplot{1,thisseason}{ensNo},var);
    rmind = CPM_fut.dataplot{1,thisseason}{ensNo}.pIMF<fpIMF;
    A(rmind) = [];
    f(ensNo,:) = ksdensity(A,yi,'function','pdf','kernel',kernelName);% 'Bandwidth',bw,
end
plot(yi,sm(nanmedian(f,1)),'b-');
hsimran = fill([yi,flip(yi)],[sm(nanmin(f,[],1)),sm(flip(nanmax(f,[],1)))],'b',...
    'LineStyle','none');

alpha(0.3)
A = getfield(RAD.dataplot{1,thisseason},var);
rmind = RAD.dataplot{1,thisseason}.pIMF<fpIMF;
A(rmind) = [];
fobs = ksdensity(A,yi,'function','pdf','kernel',kernelName);
plot(yi,sm(fobs),'k-');
box on;

hold off;

end


