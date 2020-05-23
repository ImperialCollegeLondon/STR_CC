
colsim = [1,0,0];
colobs = [0.3,0.3,0.3];
lineType = 'x';%'-';
lineTypeObs = '^';% 'o-';
linewidth = 2;
linewidthObs = 4;
Name = {'SWestuk'}%{'London','Westuk','Scotland'};%{'London','SWestuk','Westuk','Scotland'};
dataSP = 'H:\DATA_CLIMATE\UKCP18\';

REGIONS = REGIONS_info();
close all

%% Figure 3: Plot two cdfs for given aggregation area from rawIMF

for thisstation = 1:length(Name)
    
    [h,hobs] = deal([]);
    CPM = load(sprintf('%srawIMF_Area_%s_1980-2000_4Season.mat',dataSP,Name{thisstation}),'rawIMF');
    CPM_fut = load(sprintf('%srawIMF_Area_%s_2060-2080_4Season.mat',dataSP,Name{thisstation}),'rawIMF');
    RAD = load(sprintf('%sObs_rawIMF_Area_%s_4Season.mat',dataSP,Name{thisstation}),'rawIMF');
    for aggNo = [1,7,23]%[1,23,24]%1%[1,
        ha = tight_subplot(1,2,[.15 .02],[.2 .1],[.15 .05]);
        setFigureProperty('Paper');
        set(0,'defaultAxesFontSize',10,...
    'defaultTextFontSize',10);

        thisplot = 1;
        for thisseason = [2,4]
            
            axes(ha(thisplot));
            
            aggArea = CPM.rawIMF{thisseason}(1).areaAgg(aggNo);
            [~,h(thisseason),pdd_cpm] = plot_IMFCDF(CPM.rawIMF{thisseason},aggNo,...
                colsim,linewidth);
            % [~,hfut(thisseason),pdd_cpm_fut] = plot_IMFCDF(CPM_fut.rawIMF{thisseason},aggNo,...
            %    'r',linewidth);
            [~,hobs(thisseason),pdd_rad] = plot_IMFCDF(RAD.rawIMF{thisseason},aggNo,...
                colobs,linewidthObs);
            
            YL = ylim;
            XL = xlim;
            % legend(gca,[hfut(thisseason),h(thisseason),hobs(thisseason)],...
            %     [strcat(Name(thisstation),'-Future',sprintf('(%.1f%%-%.1f%%)',...
            %     100*min(pdd_cpm_fut),100*max(pdd_cpm_fut))),...
            %     strcat(Name(thisstation),'-1980-2000',sprintf('(%.1f%%-%.1f%%)',...
            %     100*min(pdd_cpm),100*max(pdd_cpm))),...
            %     strcat(Name(thisstation),'-Radar',sprintf('(%.1f%%)',...
            %     100*pdd_rad))],...
            %     'Location','NorthEast','fontsize',20);
            %         legend(gca,[h(thisseason),hobs(thisseason)],...
            %             [strcat(Name(thisstation),'-UKCP',sprintf('(%.1f%%-%.1f%%)',...
            %             100*min(pdd_cpm),100*max(pdd_cpm))),...
            %             strcat(Name(thisstation),'-Radar',sprintf('(%.1f%%)',...
            %             100*pdd_rad))],...
            %             'Location','SouthEast');
            legend(gca,[h(thisseason),hobs(thisseason)],...
                {strcat('UKCP',sprintf('(%.1f%%-%.1f%%)',...
                100*min(pdd_cpm),100*max(pdd_cpm)));...
                strcat('Radar',sprintf('(%.1f%%)',...
                100*pdd_rad))},...
                'Location','NorthEast','fontsize',8);
            % grid minor
            title([sprintf('%s(%4.1fKm^2)',getSeasonName(thisseason),aggArea)]);
            legend boxoff
            set(gca,'linewidth',2)
            
            thisplot = thisplot+1;
        end
        
        format_xylabel(ha,1,2)
        filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\SpatialAnalysis\';
        filename = [filePath,filesep,sprintf('IMF_rad_cpm_%s_%dkm2',Name{thisstation},round(aggArea))];
        savePlot(filename,'units','centimeters','XYWH',[5,0,12,7],'needreply','N');
        close all
    end
    
end

function [h1,h2,pdd] = plot_IMFCDF(rawIMF,aggNo,colv,linewidth)

h1 = 0;

[y,f_ens,x_ens,pdd] = deal([]);
x_ens = 0.01:0.5:40;%80;%
thr = 0.1;
for enNo = 1:length(rawIMF)
    thisy = rawIMF(enNo).AreaIMFs{aggNo};
    pdd(enNo) = sum(thisy<thr)./numel(thisy);
    thisy = thisy(thisy>=thr);
    y = [y,thisy];
    [f_temp,x_temp,~,~] = ecdf(thisy);
    [x_temp,IA,~] = unique(x_temp);
    f_ens(enNo,:) = interp1(x_temp,f_temp(IA),x_ens);
end

[f,x,flo,fup] = ecdf(y(:));
h2 = plot(x,1-f,'-','color',[colv,0.8],'linewidth',linewidth);
hold on;

if size(f_ens,1) > 1
    %     fill([x_ens,flip(x_ens)],[1-prctile(f_ens,5,1),flip(1-prctile(f_ens,95,1))],'r',...
    %         'LineStyle','none');hold on;
    
    plot(x_ens,1-prctile(f_ens,5,1),'--','color',[colv,0.5],'linewidth',linewidth);
    hold on;
    plot(x_ens,1-prctile(f_ens,95,1),'--','color',[colv,0.5],'linewidth',linewidth);
    hold on;
end

configPlot_imfcdf(aggNo)

end


function configPlot_imfcdf(aggNo)
set(gca,'YScale','log')

if aggNo > 10
    xlim([-0.01,15])
elseif aggNo > 5
    xlim([-0.01,25])
else
    xlim([-0.01,50])
end
XL = xlim;
grid minor
YL = ylim;
ylim([3e-4,YL(2)])

xlabel('Areal Mean Rainrate mm/h')
ylabel('Exceedence probability');

set(gca,'YTick',[1e-3,1e-2,1e-1,1e0])
set(gca,'YTickLabel',{'1e-3','1e-2','0.1','1'});

end







