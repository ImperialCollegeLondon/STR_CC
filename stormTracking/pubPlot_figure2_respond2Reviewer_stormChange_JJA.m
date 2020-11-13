% ----------------------------------------------------------------------- %
% THIS SCRIPT IS TO VISUALIZE SEVERAL STATISTICS FOR CONVECTIVE STROM DURING
% SUMMERTIME
%
% ... # some necessary description #
% All storms (convective) simulated in CPM are identified.
%
% Convective:
% Those storms having point rainfall (at least once) larger than 5mm/h
%            are identified and are regarded as convective storms.
% Several statistics for each snapshot within those identified storms were
% extracted;
% Including
%    rvol
%    rsize
%    rpmax
%    rspeed
%
% Further analyses can be done to extract those chacasteristics for each
% storm in stead for all hours.
%
% @ Yuting Chen
% Imperial College London
% yuting.chen17@imperial.ac.uk
% ----------------------------------------------------------------------- %

% Several Config
warning on

ENSEMBLENO = getEnsNos();
% 1:12
[CELLSPEED,STORMSPEED] = deal([]);
for etag = 1:12
    ensNo = ENSEMBLENO{etag};
    periodTag = 1;
    for Period = {'1980-2000','2060-2080'}
        MON = [6:8];
        Period = Period{1};% 
        regionTag = 1;
        for regionName = {'CPM_NW','CPM_NE','CPM_S'}
            regionName = regionName{1};
            
            Config = getConfig(upper(regionName),8,Period,ensNo);
            Config.Month = [MON];
            load([Config.saveIt.path,filesep,sprintf('CS_%s_%s_FIELD_%02d-%02d_%s.mat',regionName,Period,...
                                Config.Month(1),Config.Month(end),ensNo)],'RE','TE');
            
            STATS = getSTATS(RE,Config);
            cellSpeed = cell2mat(cellfun(@(x)x.cellSpeed,STATS,'UniformOutput',false));
            stormSpeed = cell2mat(cellfun(@(x)x.rspeed,STATS,'UniformOutput',false));
            cellSpeed(isnan(cellSpeed)) = [];
            stormSpeed(isnan(stormSpeed)) = [];
            %%
            histogram(cellSpeed(cellSpeed>0.1),'Normalization','pdf')
            hold on;histogram(stormSpeed,'Normalization','pdf')
            xlim([0,144])
            legend('Cells','Storm');
            xlabel('Speed [km/h]');
            ylabel('PDF')
            pause(5)
            %%
            close all
            CELLSPEED{regionTag,periodTag,etag} = cellSpeed;
            STORMSPEED{regionTag,periodTag,etag} = stormSpeed;
            regionTag = regionTag+1;
        end
        periodTag = periodTag+1;
    end
     
end

save('speed_response.mat','CELLSPEED','STORMSPEED','-v7.3')

%% obtain information required to compute CELLSPEED and STORMSPEED
load('speed_response.mat','CELLSPEED','STORMSPEED')
[N_cellSpeed,N_stormSpeed] = deal([]);
for regionID = 1:3
    for periodID = 1:2
        for ensNo = 1:12
            % subplot(3,4,ensNo)
            cellSpeed = CELLSPEED{regionID,periodID,ensNo};
            stormSpeed = STORMSPEED{regionID,periodID,ensNo};
            [N_cellSpeed(regionID,ensNo,periodID,:),EDGES] = histcounts(cellSpeed(cellSpeed>0),'BinEdges',0:2:200,'Normalization','pdf');
            hold on;[N_stormSpeed(regionID,ensNo,periodID,:),EDGES] = histcounts(stormSpeed,'BinEdges',0:2:200,'Normalization','pdf');
            % xlim([0,144])
            % legend('Cells','Storm');
            % xlabel('Speed [km/h]');
            % ylabel('PDF')
            % ylim([0,0.05])
        end
        pause(5)
        close all
    end
end
%% plot CELLSPEED histogram and STORMSPEED histogram
for periodID = 1:2
    subplot(1,2,periodID)
    X0 = squeeze(nanmean(N_stormSpeed,[1,2]));
    X1 = squeeze(nanmean(N_cellSpeed,[1,2]));
    h = bar(EDGES([1:100])+EDGES([2:101]),[X0(periodID,:);X1(periodID,:)],2);
    hold on;
    h(1).FaceColor = [0.5 0.5 0.5];h(1).FaceAlpha = 0.8;
    h(2).FaceAlpha = 0.8;
    xlim([0,144])
    legend('Storm','Cells');
    xlabel('Speed [km/h]');
    ylabel('PDF')
    ylim([0,0.05])
    axis('square');
    title(sprintf('(%s)','a'+periodID-1));
end

savePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\StormTracking';
fileName = 'stormSpeedvscellSpeed';
savePlot([savePath,filesep,fileName],'targetSize','dc','needreply','Y');% 'XYWH',[150,0,750,320],'needreply','N');

%% AUXILLARY FUNCTION

function STATS = updateSize(RE,Config,STATS,thr);
UKMap = getUKMap();
[E,N] = getEN(Config.region);
in = inpolygon(E,N,UKMap.borderE/1000,UKMap.borderN/1000);

% STATS = cellfun(@(R,i)compute4OneStorm(R,i),RE,num2cell([1:length(RE)]'),'UniformOutput',false);
rsize = cell2mat(cellfun(@(R)compute4OneStorm(R,thr),RE,'UniformOutput',false));
STATS.rsizeP999 = STATS.rsize;
STATS.rsize = rsize;

    function rsize = compute4OneStorm(R,thr)
%         m3d = double(R);
%         In = repmat(in,[1,1,size(m3d,3)]);
%         m3d(~In) = NaN;
        R = reshape(R,[],size(R,3));
        R(~in,:) = NaN;
        % R = m3d;%squeeze(m3d(:,:,rpmax>5));
        rsize = computeRSIZE(R,thr);% unit [km^2]
    end
    function rsize = computeRSIZE(R,rainThre)
        % unit: mm.*km*km/hour
        % R = squeeze(R);
        areaMat = (2.2)^2;
        rsize = nansum(R>rainThre,1);
        rsize = rsize*areaMat;
        rsize = rsize(:);
    end
end

function T = getTime(Period,Mon)
if strcmp(Period,'2007-2018')
    T = datetime(2007,1,1,0,0,0):1/24:datetime(2018,12,31,23,0,0);
    T = T(T.Month == Mon);
    T = T(:);
else
    yearRange = strcmp(Period,'2060-2080')*(2061:2080)+...
        strcmp(Period,'2020-2040')*(2021:2040)+...
        strcmp(Period,'1980-2000')*(1981:2000);
    if Mon == 12
        yearRange = yearRange-1;
    end
    [hh,dd,yy] = meshgrid(1:24,1:30,yearRange);
    T = datetime(yy,Mon,dd,hh,0,0);
    T = permute(T,[2,1,3]);
    T = T(:);
end
end

function STATS = getSTATS(RE,Config);
%    Precipitation volume
%    Size
%    Pmax
%    Speed

% trim area
UKMap = getUKMap();
[E,N] = getEN(Config.region);
in = inpolygon(E,N,UKMap.borderE/1000,UKMap.borderN/1000);

STATS = cellfun(@(R,i)compute4OneStorm(R,i),RE,num2cell([1:length(RE)]'),'UniformOutput',false);
rsizeall = cell2mat(cellfun(@(x)x.rsizeall,STATS,'UniformOutput',false));
rspeed = cell2mat(cellfun(@(x)x.rspeed,STATS,'UniformOutput',false));
cellSpeed = cell2mat(cellfun(@(x)x.cellSpeed,STATS,'UniformOutput',false));
rdur = cell2mat(cellfun(@(x)x.rdur,STATS,'UniformOutput',false));
stats = struct('rsizeall',rsizeall,'rspeed',rspeed,'rdur',rdur);
stats.cellSpeed = cellSpeed;

    function stats = compute4OneStorm(R,i)
        
        m3d = double(R);
        R = m3d;%squeeze(m3d(:,:,rpmax>5));
        rpmax = nanmax(reshape(m3d,[],size(m3d,3)),[],1);
        rpmax = rpmax(:);
        [rspeed,cellSpeed] = computeRSPEED(R);% unit [km/h]
        
        In = repmat(in,[1,1,size(m3d,3)]);
        m3d(~In) = NaN;
        R = m3d;%squeeze(m3d(:,:,rpmax>5));
        rsizeall = computeRSIZE(R,0.1);% unit [km^2]
        rdur = computeRDUR(R);% unit [h]
        stats = struct('rsizeall',rsizeall,'rspeed',rspeed,'rdur',rdur,'rpmax',rpmax);
        stats.cellSpeed = cellSpeed;
    end


    function rsize = computeRSIZE(m3d,rainThre)
        % unit: mm.*km*km/hour
        m3d = squeeze(m3d);
        areaMat = (2.2)^2;
        rsize = nansum(reshape(m3d>rainThre,[],size(m3d,3)),1);
        rsize = rsize*areaMat;
        rsize = rsize(:);
    end

    function rdur = computeRDUR(m3d)
        rdur = repmat(size(m3d,3),[1,size(m3d,3)]);
        rdur = rdur(:);
    end
    function [rspeed,cellSpeed] = computeRSPEED(m3d)
        
        % imreso = 0.2;% aggregate to 11km.
        
        rspeed = [];cellSpeed = [];
        zeroVal = 0;
        Rtemp = permute(m3d,[3,1,2]);% make time steps in the first dimension;
        rpmax = squeeze(nanmax(m3d,[],[1,2]));
        
        Rtemp = func_onlyMCS(Rtemp,5,zeroVal);
        
        opticFlow = opticalFlowLK('NoiseThreshold',0.0005);
        
        R0 = squeeze(Rtemp(1,:,:));%imresize(,imreso,'box');
        R0 = conv2(R0,ones(50)/50/50,'same');
        R0 = func_R2dBZ(R0,'UKMO',zeroVal);
        for sni = 1:size(Rtemp,1)
            
            R1 = squeeze(Rtemp(sni,:,:));
            R1 = conv2(R1,ones(50)/50/50,'same');
            R1 = func_R2dBZ(R1,'UKMO',zeroVal);
            
            % mx = nanmax([R1(R1>zeroVal);R2(R2>zeroVal);zeroVal]);
            % mn = nanmin([R1(R1>zeroVal);R2(R2>zeroVal)]);
            % if isempty(mn)
            %     mn = 0;
            % end
            % R1 = (R1-mn)/(mx-mn);R1(R1<zeroVal) = zeroVal;
            % R2 = (R2-mn)/(mx-mn);R2(R2<zeroVal) = zeroVal;
            
            flow = estimateFlow(opticFlow,R1);
            
            V = 2.2 * nanmean(flow.Magnitude(R1~=zeroVal | R0~=zeroVal));
            rspeed(sni) = V;
            stats = regionprops(logical(func_R2dBZ(squeeze(Rtemp(sni,:,:)),'UKMO',zeroVal)),'PixelIdxList');
            cellSpeed{sni} = cell2mat(arrayfun(@(x)2.2* nanmean(flow.Magnitude(x.PixelIdxList)),stats,...
                'UniformOutput',false));
            if (~any(R0(:)~=zeroVal) & any(R1(:)~=zeroVal))|...
                    (~any(R1(:)~=zeroVal) & any(R0(:)~=zeroVal))%#ok<OR2,AND2>
                % In the case that: having peak(>5mm/h) pixels in image 0/1 but no peak
                % in image 1/0, the speed is unknown (because the things
                % happen within one hour which is not able to be captured
                % in current output.
                rspeed(sni) = NaN;
                cellSpeed{sni} = [];
            end
            hold off;
            R0 = R1;
        end
        
        rspeed = rspeed(rpmax>10);
        cellSpeed = cellSpeed(rpmax>10);
        
        rspeed = rspeed(:);
        cellSpeed = cell2mat(cellSpeed');
        
        function Rtemp = func_onlyMCS(Rtemp,thres,zeroVal);
            Rtemp(Rtemp<thres) = zeroVal;
        end
        function Rtemp = func_R2dBZ(Rtemp,source,zeroVal)
            if strcmp(source,'UKMO')
                Rtemp(Rtemp~=zeroVal) = (log10(200)+1.6*log10(Rtemp(Rtemp~=zeroVal)))*10;
            end
        end
    end
end

