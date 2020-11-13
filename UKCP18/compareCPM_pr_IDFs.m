%% read CPM_pr point
% This script is to compare the IDFs computed from the CPM2.2 and MIDAS rain guages.
% it mainly includes the following sections:
% 1. getRG information
% 2. read the .nc files to extract those points
% 3. compute IDFs for those locations having data from both datasets.
% 4. visualize IDFs in several different ways

%% get all available hourly rain gauges operating earlier than 1980 and later than 2000;
function compareCPM_pr_IDFs()
T = getRGInfo();
save('D:\UKCP18\RG\midas_rd.mat','T');

% %%
% ncFileName = ['K:/UkCp18/01/pr_rcp85_land-cpm_uk_2.2km_01_1hr_19910801-19910830.nc'];
% LAT=ncread(ncFileName,'latitude');LON=ncread(ncFileName,'longitude');
% [E, N] = ll2os(LAT, LON);[e0,n0] = ll2os(T.Latitude,T.Longitude);
% E = E/1000; N = N/1000;
% e0 = e0/1000; n0 = n0/1000;
% [RGi,RGj] = arrayfun(@(e0,n0)getRegionIJ(E,N,e0,n0),e0,n0);
% 
% ENSEMBLENO = getEnsNos();
% MON = 6:8;
% data = struct('Years',[1980,2000],'fileGetPath','K:/UkCp18/',...
%     'savePath','D:/UKCP18/RG',...
%     'RGi',RGi,'RGj',RGj);
% options = 'Points';
% [RainEnsembles] = readCPM_pr_annualMax(ENSEMBLENO,MON,data,options);
% 
% %% Compute IDF of CPM for each site
% duration = hours([1:24]);
% returnT = [2:20];
% IDFEnsembles = [];
% for ensNo = 1:12
%     intensity = [];
%     returnT = reshape(returnT,[],1);
%     A = load(['D:\UKCP18\RG\ensembles',num2str(ensNo),'_RG.mat']);
%     rain = A.RainEnsembles(1,2:21);clear A
%     rain = cell2mat(cellfun(@(a)reshape(a,[1,size(a)]),rain',...
%         'UniformOutput',false));
%     index = 1;
%     for dur = duration
%         tic
%         rain0 = imresize3(rain,'scale',[1,1,1/hours(dur)],'method','box');
%         am = getAM(rain0);
%         for sitei = 1:size(am,2)
%             intensity(index,:,sitei) = getReturnLevels(am(:,sitei),returnT);
%         end
%         index = index+1;
%         toc
%     end
%     IDFEnsembles{ensNo} = intensity;
% end
% save(['D:\UKCP18\RG\ensembles',num2str(ensNo),'_RG_IDF.mat'],'IDFEnsembles',...
%     'duration','returnT');
% 
% 
% %% Compute IDF of MIDAS for each site
% furation = hours([1:24]);
% returnT = [2:20];
% 
% intensity = [];
% returnT = reshape(returnT,[],1);
% [rain,time] = getMIDAS(T);
% 
% %%
% index = 1;
% for sitei = 1:size(am,2)
%     thisRain = rain{sitei};thisTime = datetime(datevec(time{sitei}));
%     if ~isempty(thisRain)
%         thisRain = thisRain(ismember(thisTime.Month,[6,7,8]));
%         thisRain = reshape(thisRain,nansum(eomday(2020,[6,7,8]))*24,[]);
%         thisRain = permute(thisRain,[2,1]);
%         rain{sitei} = thisRain;
%     else
%         rain{sitei} = NaN(21,nansum(eomday(2020,[6,7,8]))*24);
%     end
% end
% rain = cell2mat(cellfun(@(x)reshape(x,[1,size(x)]),rain','UniformOutput',false));
% rain = permute(rain,[2,1,3]);
% 
% for dur = duration
%     tic
%     rain0 = imresize3(rain,'scale',[1,1,1/hours(dur)],'method','box');
%     am = getAM(rain0);
%     for sitei = 1:size(am,2)
%         if isnan(am(1,sitei))
%         else
%             intensity(index,:,sitei) = getReturnLevels(am(:,sitei),returnT);
%         end
%     end
%     index = index+1;
%     toc
% end
% IDFEnsembles{1} = intensity;
% save(['D:\UKCP18\RG\midas_RG_IDF.mat'],'IDFEnsembles',...
%     'duration','returnT');

load(['D:\UKCP18\RG\midas_RG_IDF.mat'],'IDFEnsembles',...
    'duration','returnT');
CPM = load(['D:\UKCP18\RG\ensembles12_RG_IDF.mat'],'IDFEnsembles',...
    'duration','returnT');
% Plot IDF
[Tind] = plotIDFs_eachRainGauges();
%%
% Plot RG sites over the whole UK (~43 sites)
[E,N] = plotRG_location();
pause(3);close all
% plot IDF for each sub-region separately
%
plotIDFs_eachRegion();

% plot IDF for all RGs in one region

plotIDFs_RGsInOneRegion();

% AUXILLARY FUNCTION
    function plotIDFs_eachRegion()
        seeDR = [1,3,6,12,24];
        pl = 1;
        ha = tight_subplot(3,3,[.03 .03],[.10 .1],[.10 .05]);
        setFigureProperty('Paper')
        for regionName = {'NWUK','NEUK','SUK'}
            regionName = regionName{1};
            inRG = getInThisRegion(regionName,E(Tind)/1000,N(Tind)/1000);
            for rti = [1,4,9]
                axes(ha(pl))
                x = hours(CPM.duration);
                idf = squeeze(IDFEnsembles(:,rti,Tind(inRG)));
                x = x(seeDR);idf = idf(seeDR,:);
                X = [x,flip(x)];
                Y = [prctile(idf,0,2);flip(prctile(idf,100,2))];% Y = exp(smooth(log(Y),'sgolay'));
                hobsR = fill(X,Y,'r','LineStyle','none','facealpha',0.3);hold on
                idf = nanmean(idf,2);
                hobs = plot(x,idf,'-','linewidth',1,'color',[1 0 0 0.4],'linewidth',2);hold on
                
                IDF = cell2mat(cellfun(@(x)x(:,rti,Tind(inRG)),CPM.IDFEnsembles,'UniformOutput',false));
                IDF = IDF(seeDR,:,:);
                IDF = squeeze(nanmedian(IDF,2));
                Y = [prctile(IDF,0,2);flip(prctile(IDF,100,2))];
                hcpmRs = fill(X,Y,'k','LineStyle','none','facealpha',0.3);hold on
                IDF = cell2mat(cellfun(@(x)x(:,rti,Tind(inRG)),CPM.IDFEnsembles,'UniformOutput',false));
                IDF = IDF(seeDR,:,:);
                IDF = reshape(IDF,[size(IDF,1)],[]);
                Y = [prctile(IDF,100*0.001/12,2);flip(prctile(IDF,100*11.999/12,2))];
                hcpmR12 = fill(X,Y,'k','LineStyle','none','facealpha',0.2);hold on
                IDF = nanmedian(IDF,2);
                hcpm = plot(x,IDF,'-','linewidth',1,'color',[0 0 0 0.4],'linewidth',2);hold on
                hold on;
                set(gca,'XScale','log')
                set(gca,'YScale','log')
                xlim([1,24]);ylim([1,40]);ax = gca;ax.XTick = [1,2,6,12,24];ax.YTick = [2,5,10,20];
                
                legend([hcpm,hobs,hcpmRs,hcpmR12,hobsR],'CPM2.2','Obs','CPM-RangeSites','CPM-RangeEns','Obs-Range','fontsize',10)
                xlabel('duration [h]','fontsize',10);
                ylabel('intensity [mm/h]','fontsize',10);
                str = sprintf('RT:%dyr-%s',rti+1,regionName);
                text(23.5,35,str,'horizontalAlignment','right','VerticalAlignment','top',...
                    'fontsize',8);
                grid minor
                % axis('square')
                savePlot(['temp'],'targetSize','1c','needreply','N','onlyPng',true);
                pl = pl+1;
            end
        end
        format_xylabel(ha,3,3)
    end

    function plotIDFs_RGsInOneRegion()
        seeDR = [1,3,6,12,24];
        for regionName = {'SUK'};%{'NWUK','NEUK','SUK'}
            regionName = regionName{1};
            inRG = getInThisRegion(regionName,E(Tind)/1000,N(Tind)/1000);
            figure;
            setFigureProperty('Paper')
            rti = 9;
            pl = 1;
            ha = tight_subplot(5,5,[.02 .02],[.05 .05],[.10 .05]);
            for ti = Tind(inRG)
                axes(ha(pl))
                x = hours(CPM.duration);
                idf = squeeze(IDFEnsembles(:,rti,ti));
                x = x(seeDR);idf = idf(seeDR);
                X = [x,flip(x)];
                Y = [prctile(idf,0,2);flip(prctile(idf,100,2))];
                hobsR = plot(X,Y,'r');hold on
                
                
                IDF = cell2mat(cellfun(@(x)x(:,rti,ti),CPM.IDFEnsembles,'UniformOutput',false));
                IDF = IDF(seeDR,:);
                Y = [prctile(IDF,100*1.1/12,2);flip(prctile(IDF,100*10.9/12,2))];
                hcpmR12 = fill(X,Y,'k','LineStyle','none','facealpha',0.2);hold on
                IDF = nanmedian(IDF,2);
                hcpm = plot(x,IDF,'-','linewidth',1,'color',[0 0 0 0.4],'linewidth',2);hold on
                hold on;
                set(gca,'XScale','log')
                set(gca,'YScale','log')
                xlim([1,25]);ylim([1,40]);ax = gca;ax.XTick = [1,2,6,12,24];ax.YTick = [2,5,10,20];
                
                % legend([hcpm,hobs,hcpmRs,hcpmR12,hobsR],'CPM2.2','Obs','CPM-RangeSites','CPM-RangeEns','Obs-Range','fontsize',10)
                xlabel('duration [h]','fontsize',10);
                ylabel('intensity [mm/h]','fontsize',10);
                str = sprintf('%s RGs:%d',regionName,nansum(inRG));
                % text(23.5,29,str,'horizontalAlignment','right','VerticalAlignment','top','fontsize',10);
                grid minor
                axis('square')
                pl = pl+1;
            end
            format_xylabel(ha,5,5)
        end
    end

    function [E,N] = plotRG_location()
        figure
        [E,N] = ll2os(T.Latitude,T.Longitude);
        UKMap = getUKMap();
        plot(UKMap.borderE/1000,UKMap.borderN/1000,'k-','color',[0.5 0.5 0.5 0.8]);hold on
        plot(E(Tind)/1000,N(Tind)/1000,'ko','markerfacecolor','w');
        xlim([0,800]);ylim([0,1300])
        cmap = cptcmap('flood_blue','mapping','scaled','ncol',5,'flip',true);
        
        angle = 15;
        
        regionName = 'SUK';
        leng = 88;
        ul_x = 50;ul_y = 250;
        ur_x = 660;ur_y = ul_y+tan(angle/180*pi)*(ur_x-ul_x);
        X = [ul_x,ul_x+leng,ur_x+leng,ur_x,ul_x];
        Y = [ul_y,ul_y-(leng/tan(angle/180*pi)),ur_y-(leng/tan(angle/180*pi)),ur_y,ul_y];
        oneRegion(X,Y,regionName,cmap(1,:))
        
        regionName = 'NWUK';
        ll_x = X(1);ll_y = Y(1);
        lr_x = 400;lr_y = ll_y+tan(angle/180*pi)*(lr_x-ll_x);
        leng = 200;
        X = [ll_x-leng,ll_x,lr_x,lr_x-leng,ll_x-leng];
        Y = [ll_y+(leng/tan(angle/180*pi)),ll_y,lr_y,lr_y+(leng/tan(angle/180*pi)),ll_y+(leng/tan(angle/180*pi))];
        oneRegion(X,Y,regionName,cmap(2,:))
        
        regionName = 'NEUK';
        ll_x = X(3);ll_y = Y(3);
        lr_x = 600;lr_y = ll_y+tan(angle/180*pi)*(lr_x-ll_x);
        leng = 200;
        X = [ll_x-leng,ll_x,lr_x,lr_x-leng,ll_x-leng];
        Y = [ll_y+(leng/tan(angle/180*pi)),ll_y,lr_y,lr_y+(leng/tan(angle/180*pi)),ll_y+(leng/tan(angle/180*pi))];
        oneRegion(X,Y,regionName,cmap(3,:))
        
        hold off
        axis equal
        xlim([0,800]);
        
        ax = gca;
        ax.XTick = [];
        ax.YTick = [];
        ax.LineWidth = 0.5;
        box on
        
        % savePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\';
        % savePlot('temp','targetSize','1c','needreply','Y','onlyPng',true);
    end

    function T = getRGInfo()
        file = 'H:\DATA_CEDA\src_ids.xlsx';
        T = readtable(file);
        
        tag = cellfun(@(strStart,strEnd)str2num(strStart(end-3:end))<1980 & ...
            (strcmp(strEnd,'Current')||str2num(strEnd(end-3:end))>2000),T.StationStartDate,T.StationEndDate);
        T = T(tag,:);
    end

    function [Tind] = plotIDFs_eachRainGauges();
        rti = 2;
        Tind = [];
        itag = 1;
        ha = tight_subplot(5,9,[.03 .01],[.05 .05],[.10 .05]);
        x = hours(CPM.duration);
        for sitei = 1:size(IDFEnsembles,3)
            idf = IDFEnsembles(:,rti,sitei);
            if idf(:,1)>40
                sitei;
            end
            if idf(:,1)~=0 & idf(:,1)<40 %#ok<AND2> %
                Tind(itag) = sitei;
                axes(ha(itag))
                IDF = cell2mat(cellfun(@(x)x(:,rti,sitei),CPM.IDFEnsembles,'UniformOutput',false));
                X = [x,flip(x)];
                Y = [prctile(IDF,100*0.1/12,2);flip(prctile(IDF,100*11.9/12,2))];
                Y = exp(smooth(log(Y),'sgolay'));X = exp(smooth(log(X),'sgolay'));
                hcpm = fill(X,Y,[0.5 0.5 0.5],'LineStyle','none');
                hold on;
                hobs = plot(x,idf,'r-','linewidth',1);
                hold on;
                set(gca,'XScale','log')
                set(gca,'YScale','log')
                xlim([1,24]);ylim([2,30]);ax = gca;ax.XTick = [1,2,6,12,24];ax.YTick = [2,5,10,20];
                itag = itag+1;
            end
        end
        legend([hcpm,hobs],'cpm2.2','obs')
    end

    function in = getInThisRegion(thisRegionName,E,N)
        angle = 15;
        
        regionName = 'SUK';
        leng = 88;
        ul_x = 50;ul_y = 250;
        ur_x = 660;ur_y = ul_y+tan(angle/180*pi)*(ur_x-ul_x);
        regionX = [ul_x,ul_x+leng,ur_x+leng,ur_x,ul_x];
        regionY = [ul_y,ul_y-(leng/tan(angle/180*pi)),ur_y-(leng/tan(angle/180*pi)),ur_y,ul_y];
        in = inpolygon(E,N,regionX,regionY);
        if strcmp(regionName,thisRegionName)
            return
        end
        regionName = 'NWUK';
        ll_x = regionX(1);ll_y = regionY(1);
        lr_x = 400;lr_y = ll_y+tan(angle/180*pi)*(lr_x-ll_x);
        leng = 200;
        regionX = [ll_x-leng,ll_x,lr_x,lr_x-leng,ll_x-leng];
        regionY = [ll_y+(leng/tan(angle/180*pi)),ll_y,lr_y,lr_y+(leng/tan(angle/180*pi)),ll_y+(leng/tan(angle/180*pi))];
        in = inpolygon(E,N,regionX,regionY);
        if strcmp(regionName,thisRegionName)
            return
        end
        regionName = 'NEUK';
        ll_x = regionX(3);ll_y = regionY(3);
        lr_x = 600;lr_y = ll_y+tan(angle/180*pi)*(lr_x-ll_x);
        leng = 200;
        regionX = [ll_x-leng,ll_x,lr_x,lr_x-leng,ll_x-leng];
        regionY = [ll_y+(leng/tan(angle/180*pi)),ll_y,lr_y,lr_y+(leng/tan(angle/180*pi)),ll_y+(leng/tan(angle/180*pi))];
        in = inpolygon(E,N,regionX,regionY);
        if strcmp(regionName,thisRegionName)
            return
        end
    end

    function oneRegion(X,Y,regionName,color)
        line(X,Y,'color',color);
        hold on;
        [x0,y0] = centroid(polyshape({X},{Y}));
        text(x0,y0,regionName,...
            'fontsize',12,'horizontalalignment','center','fontweight','bold','background',[0.9 0.9 0.9 0.6]);
    end

    function [rain,time] = getMIDAS(T)
        % rain = [];
        % time = [];
        % for sitei = 1:size(T,1)
        %     filename = ['H:\DATA_CEDA\DataInDifferentStation_New discardLongObservationCount\station',num2str(T.src_id(sitei)),'.txt'];
        %     fid = fopen(filename);
        %     if fid ~= -1
        %         A = load(filename);
        %         time{sitei} = datetime(datevec(A(:,1)));
        %         if time{sitei}(1,:).Year<1981 & time{sitei}(end,:).Year>2000
        %         rain{sitei} = A(:,3);
        %         fprintf('site %04d: ',sitei);
        %         [time{sitei},rain{sitei}] = interpRG(time{sitei},rain{sitei});
        %         end
        %     else
        %         rain{sitei} = [];
        %         time{sitei} = [];
        %     end
        % end
        % save('D:\UKCP18\RG\midas_timeseries.mat','rain','time','T')
        load('D:\UKCP18\RG\midas_timeseries.mat','rain','time','T')
    end

    function [date,depth] = interpRG(time,rain)
        % [depth,date] = interpRain(oridate,oridepth,minDate,maxDate,reso)
        reso = 60;
        ii = find(time.Year == 1980,1);
        time = time(ii:end);rain = rain(ii:end);
        oridepth = rain;
        oridate = datenum(time);
        [oridate,u] = unique(oridate);oridepth = oridepth(u);
        minDate = datenum(datetime(1980,1,1));
        maxDate = datenum(datetime(2000,12,31));
        
        if isempty(oridepth) || nanmax(time.Year)<2000 || length(oridepth)<1000
            [date,depth] = deal([]);
            return
        end
        
        maxdiff = max(diff(oridate));
        if maxdiff>10
            ind2 = find(diff(oridate) == maxdiff);
            % datetime(datevec(oridate(ind2))) %#ok<FNDSB>
            if maxdiff>300
                [date,depth] = deal([]);
                return
            end
            fprintf('Pay Attention! poor quality data: max diff:%dday\n',maxdiff);
        else
            fprintf('\n');
        end
        
        %create series for calibration process
        depth = interp1(oridate,oridepth,[minDate:1/24/(60/reso):maxDate]');
        date = [minDate:1/24/(60/reso):maxDate]';
        threshold = 0;%0.2;% daily threshold;
        depth(depth<=threshold/24/(60/reso)) = 0;
    end

    function [intensity,returnT,duration] = getIDF(oneSite,dt)
        % dt: <duration>
        % oneSite: <double> or <RainfallDataClass>
        
        if strcmp(class(oneSite),'RainfallDataClass')
            % time = oneSite.Time;
            % oneSite = originalData(oneSite);
        end
        
        returnT = [5;20;100];
        duration = [minutes([5,30]),hours([1:24])];
        
        returnT = reshape(returnT,[],1);
        index = 1;
        intensity = [];
        
        for dur = duration
            rainfallTS = aggregate(oneSite,1,dur/dt,'mm/h');
            time = rainfallTS.Time;
            % rainfallTS = squeeze(originalData(rainfallTS));
            rainfallTS = aggregate(squeeze(originalData(oneSite)),dur/dt,'mean');% use this one to prevent issue of imresize3 in aggregation
            am = getAM(rainfallTS,time);
            intensity(:,index) = getReturnLevels(am,returnT);
            index = index+1;
        end
        
    end

    function am = getAM(rain)
        am = nanmax(rain,[],3);
    end

    function intensity = getReturnLevels(am,rt)
        options = statset('gevfit');
        options.MaxIter = 500;
        P = 1-(1./rt);
        warning off
        [parmhat,parmci] = gevfit(am,0.05,options);
        intensity = gevinv(P,parmhat(1),parmhat(2),parmhat(3));
        [warnMsg,warnId] = lastwarn;
        if ~isempty(warnMsg)
            warnMsg = '';
            [parmhat,parmci] = evfit(am);%
            intensity = evinv(P,parmhat(1),parmhat(2));
        end
    end

    function [RainEnsembles] = readCPM_pr_annualMax(ENSEMBLENO,MON,data,options)
        % READCPM_PR() gives .....
        %
        % Input:region
        %       ENSEMBLENO: (1,>=1) <cell>
        %       MON:
        %       data:
        %       options:
        % Output:RainEnsembles
        
        
        % Configuration
        filePath = ['K:/UkCp18/',ENSEMBLENO{1},'/pr_rcp85_land-cpm_uk_2.2km_'];
        ncFileName = [filePath,ENSEMBLENO{1},'_1hr_19910801-19910830.nc'];
        LAT=ncread(ncFileName,'latitude');
        LON=ncread(ncFileName,'longitude');
        [E, N] = ll2os(LAT, LON);
        E = E/1000;
        N = N/1000;
        
        mkdir(data.savePath);
        fprintf(sprintf('%04d-%04d data will be extracted',data.Years(1),data.Years(end)));
        
        warning off
        
        
        for M=2:length(ENSEMBLENO)%%%%%%%%%%%%%%%%%%%%
            
            saveDir = [data.savePath,'/ensembles',num2str(M),'_RG'];
            
            RainEnsembles = [];
            filePath = [data.fileGetPath,ENSEMBLENO{M},'/'];
            
            L = 1;
            
            for year = data.Years(1):data.Years(end)
                thisYear = [];
                tic
                for mon = MON
                    thisMonth = [];
                    if ~((year==data.Years(1) && mon ~=12) || (year == data.Years(end) && mon == 12))
                        
                        fileName = sprintf('pr_rcp85_land-cpm_uk_2.2km_%s_1hr_%04d%02d01-%04d%02d30.nc',...
                            ENSEMBLENO{M},year,mon,year,mon);
                        
                        thisNCFile = dir([filePath,fileName]);
                        
                        if ~isempty(thisNCFile)
                            thisNCFile = fullfile({thisNCFile.folder},{thisNCFile.name});
                            thisNCFile = thisNCFile{1};
                            A = ncinfo(thisNCFile);
                            Rain_this = squeeze(ncread(thisNCFile,'pr',[1,1,1,1],[inf,inf,inf,inf]));
                            if strcmpi(options,'Points')
                                IND = sub2ind(size(Rain_this,[1,2]),data.RGi,data.RGj);
                                Rain_this = reshape(Rain_this,[prod(size(Rain_this,[1,2])),size(Rain_this,3)]);
                                thisMonth = Rain_this(IND,:);
                            end
                        end
                        thisYear = cat(2,thisYear,thisMonth);
                    end
                end
                RainEnsembles{L}=thisYear;
                L = L+1;
                toc
            end
            fprintf('Ensembles %02d\n',M);
            save([saveDir,'.mat'],'RainEnsembles');
        end
        
    end
end
