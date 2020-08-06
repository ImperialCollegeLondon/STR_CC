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
% :12
for ensNo = ENSEMBLENO([1:12]);%{'RAD'}%
    ensNo = ensNo{1};% 
    for Period = {'1980-2000','2060-2080'}%{'2007-2018'}%'2020-2040',
        MON = [6:8];
        Period = Period{1};% 'CPM_NW',
        for regionName = {'CPM_NE','CPM_S'}%{'EUK','SCO','WAL'}
            regionName = regionName{1};
            
%             RE = [];
%             %%
%             % Load Raw Data
%             [RainEnsembles,T,Config] = getData(regionName,MON,Period,ensNo);
%             
%             % GET JJA-CS
%             [pointer] = find_CS(RainEnsembles{1},Config);
%             plot(logical(pointer))
%             
%             % GET RFIELD for JJA-CS
%             % pointer(pointer==0) = NaN;
%             [RE,TE] = get_CS(RainEnsembles{1},T,pointer,Config);
%             save([Config.saveIt.path,filesep,sprintf('CS_%s_%s_FIELD_%02d-%02d_%s.mat',regionName,Period,...
%                 Config.Month(1),Config.Month(end),ensNo)],'RE','TE','Config')%,'-v7.3'
%             [warnMsg, warnId] = lastwarn;
%             if ~isempty(warnMsg)
%                 save([Config.saveIt.path,filesep,sprintf('CS_%s_%s_FIELD_%02d-%02d_%s.mat',regionName,Period,...
%                     Config.Month(1),Config.Month(end),ensNo)],'RE','TE','Config','-v7.3')
%             end
%             warning('')
            %%
            Config = getConfig(upper(regionName),8,Period,ensNo);
            Config.Month = [MON];
            load([Config.saveIt.path,filesep,sprintf('CS_%s_%s_FIELD_%02d-%02d_%s.mat',regionName,Period,...
                                Config.Month(1),Config.Month(end),ensNo)],'RE','TE');
            
            
            
            % Get STATS table for JJA-CS
            % including:
            %    Precipitation volume
            %    Size
            %    Pmax
            %    Speed
            % # Notice # here only the region of each study area is used, which means
            % the area beyond the study area was not considered to contribute the
            % storm.
            STATS = getSTATS(RE,Config);
            save([Config.saveIt.path,filesep,sprintf('CS_%s_%s_STATS_%02d-%02d_%s.mat',regionName,Period,...
                Config.Month(1),Config.Month(end),ensNo)],'STATS','Config')
            histogram(STATS.rspeed(STATS.rpmax>=10));
            
        end
        
    end
    
    % get STATS trajectory
    % R = double(RainEnsembles{ensNo}(:,:,1:end))/IntFac;
    % R = permute(R,[3,1,2]);
    
    
end

% AUXILLARY FUNCTION

function [RainEnsembles,T,Config] = getData(regionName,MON,Period,ENSNO)

RainEnsembles = [];
T = [];

for mon = MON(:)'
    Config = getConfig(upper(regionName),mon,Period,ENSNO);
    datafile = Config.data.file;
    A = load(datafile);RainEnsembles_temp = [];
    if isfield(A,'RainEnsembles')
        RainEnsembles_temp = A.RainEnsembles;
    else
        RainEnsembles_temp{1} = A.Rain;
    end
    T_temp = getTime(Period,mon);
    for ensNo = 1:length(RainEnsembles_temp)
        % # only valid for UKCP hourly projections #
        shapeSize = [size(RainEnsembles_temp{ensNo},[1,2]),30*24,size(RainEnsembles_temp{ensNo},3)/30/24];
        if iscell(RainEnsembles) && ~isempty(RainEnsembles{ensNo})
            RainEnsembles{ensNo} = cat(3,RainEnsembles{ensNo},...
                reshape(RainEnsembles_temp{ensNo},shapeSize));
        else
            RainEnsembles{ensNo} = reshape(RainEnsembles_temp{ensNo},shapeSize);
        end
    end
    
    T_temp = reshape(T_temp,30*24,numel(T_temp)/30/24);
    fprintf('Output is int R with scaleF\n');
    if isempty(T)
        T = T_temp;
    else
        T = [T;T_temp];
    end
end

Config.Month = MON;
RainEnsembles = cellfun(@(x)reshape(x,[size(x,[1,2]),prod(size(x,[3,4]))]),...
    RainEnsembles,'UniformOutput', false);
T = T(:);

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

function [pointer] = find_CS(R,Config)
% # output pointer (indics of R) of CSs #
%
% Currently, Convective storms were identified using following processes.
% Please notice that the first step (event seperation) is not very
% necessary for current CSs I identified.
% But it might be needed if some track algorithms are plugged in. So I keep
% this function in this way.
%

% seperate & find all events
[eventTag] = seperateEvent(R,Config.method.eventSeperation,Config.data.IntFac);

% filter out non-convective events
[eventTag] = onlyMCS(R,eventTag,Config.method.MCScrit,Config.data.IntFac);

% get pointer for those events meeting our requirement
pointer = eventTag;

end

function [RE,TE] = get_CS(R,T,pointer,Config)

[RE,TE] = deal(cell(numel(unique(pointer))-1,1));
tag = 1;
for evi = reshape(unique(pointer),1,[])
    if evi~=0
        TE{tag} = T(pointer==evi);
        RE{tag} = double(squeeze(R(:,:,pointer==evi)))./Config.data.IntFac;
        RE{tag} = single(RE{tag});
        tag = tag+1;
    end
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
rvol = cell2mat(cellfun(@(x)x.rvol,STATS,'UniformOutput',false));
rrmi = cell2mat(cellfun(@(x)x.rrmi,STATS,'UniformOutput',false));
rsize = cell2mat(cellfun(@(x)x.rsize,STATS,'UniformOutput',false));
rsizeall = cell2mat(cellfun(@(x)x.rsizeall,STATS,'UniformOutput',false));
rpmax = cell2mat(cellfun(@(x)x.rpmax,STATS,'UniformOutput',false));
rspeed = cell2mat(cellfun(@(x)x.rspeed,STATS,'UniformOutput',false));
rdur = cell2mat(cellfun(@(x)x.rdur,STATS,'UniformOutput',false));
evi = cell2mat(cellfun(@(x)x.evi,STATS,'UniformOutput',false));
STATS = table(rvol,rrmi,rsize,rsizeall,rpmax,rspeed,rdur,evi);
STATS.mon = STATS.evi*0+Config.Month;

    function stats = compute4OneStorm(R,i)
        
        m3d = double(R);
        R = m3d;%squeeze(m3d(:,:,rpmax>5));
        rspeed = computeRSPEED(R);% unit [km/h]
        
        In = repmat(in,[1,1,size(m3d,3)]);
        m3d(~In) = NaN;
        R = m3d;%squeeze(m3d(:,:,rpmax>5));
        rpmax = nanmax(reshape(m3d,[],size(m3d,3)),[],1);
        rpmax = rpmax(:);
        rvol = computeRVOL(R); % unit [m^3/s]
        rrmi = computeRMI(R); % unit [mm/h]
        rsize = computeRSIZE(R,5);% unit [km^2]
        rsizeall = computeRSIZE(R,0.1);% unit [km^2]
        rpmax = computeRPMAX(R);% unit [mm/h]
        
        rdur = computeRDUR(R);% unit [h]
        stats = table(rvol,rrmi,rsize,rsizeall,rpmax,rspeed,rdur);
        stats.evi = repmat(i,size(rvol));
    end


    function rvol = computeRVOL(m3d)
        % unit: m^3/s
        m3d = squeeze(m3d);
        areaMat = (2.2)^2;
        rvol = nansum(reshape(m3d.*areaMat,[],size(m3d,3)),1);
        rvol = rvol/1000*1000^2/3600;
        rvol = rvol(:);
    end
    function rmi = computeRMI(m3d)
        % unit: mm/h
        m3d = squeeze(m3d);
        m2d = reshape(m3d,[],size(m3d,3));
        rmi = nansum(m2d,1)./nansum(m2d>0.1,1);
        rmi = rmi(:);
    end
    function rsize = computeRSIZE(m3d,rainThre)
        % unit: mm.*km*km/hour
        m3d = squeeze(m3d);
        areaMat = (2.2)^2;
        rsize = nansum(reshape(m3d>rainThre,[],size(m3d,3)),1);
        rsize = rsize*areaMat;
        rsize = rsize(:);
    end

    function rpmax = computeRPMAX(m3d)
        rpmax = nanmax(reshape(m3d,[],size(m3d,3)),[],1);
        rpmax = rpmax(:);
    end
    function rdur = computeRDUR(m3d)
        rdur = repmat(size(m3d,3),[1,size(m3d,3)]);
        rdur = rdur(:);
    end
    function rspeed = computeRSPEED(m3d)
        
        % imreso = 0.2;% aggregate to 11km.
        
        rspeed = [];
        zeroVal = 0;
        Rtemp = permute(m3d,[3,1,2]);% make time steps in the first dimension;
        
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
            % fprintf('%2.2d\n',V);
            
            if (~any(R0(:)~=zeroVal) & any(R1(:)~=zeroVal))|...
                    (~any(R1(:)~=zeroVal) & any(R0(:)~=zeroVal))%#ok<OR2,AND2>
                % In the case that: having peak(>5mm/h) pixels in image 0/1 but no peak
                % in image 1/0, the speed is unknown (because the things
                % happen within one hour which is not able to be captured
                % in current output.
                rspeed(sni) = NaN;
            end
            hold off;
            R0 = R1;
        end
        
        rspeed = rspeed(:);
        % plot(rspeed);
        % drawnow
        % pause(0.2);
        
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



%%
function [eventTag] = onlyMCS(R,eventTag,method,scaleF);
switch(method)
    case 'PIMF>5'
        PIMF = double(nanmax(reshape(R,[],size(R,3)),[],1))/scaleF;
        T = table(PIMF(:),eventTag(:),...
            'VariableNames',{'PIMF','EvNo'});
        A = grpstats(T,'EvNo',{'max','numel'});
        eventTag = T.EvNo;
        for exEv = reshape(A.EvNo(A.max_PIMF<5 | A.numel_PIMF<0),1,[])
            eventTag(eventTag == exEv) = 0;
        end
        % #TO DO# maybe threshold for duration will be required at some
        % point.
    case 'All'
        PIMF = double(nanmax(reshape(R,[],size(R,3)),[],1))/scaleF;
        T = table(PIMF(:),eventTag(:),...
            'VariableNames',{'PIMF','EvNo'});
        A = grpstats(T,'EvNo',{'max','numel'});
        eventTag = T.EvNo;
end
end
function [eventTag] = seperateEvent(R,method,scaleF)

arguments
    R (:,:,:)
    method (1,:) char = 'WAR-based'
    scaleF (1,1) double = 32
end

switch(method)
    case 'WAR-based'
        [WAR] = compute_WetArea(R,scaleF,0.1);
        [lengthi,stormsi] = seperateRadarEvents( WAR, 0.02, 2, 0, 60);
        EvNo = [zeros(1,stormsi(1)-1),cell2mat(arrayfun(@(si,nsi,li,evi)...
            [evi*ones(1,li),0*ones(1,nsi-si-li)],stormsi,...
            [stormsi(2:end),length(WAR)+1],lengthi,1:length(stormsi),'UniformOutput', false))];
    otherwise
        error('CHECK');
end

eventTag = EvNo;

    function [WAR] = compute_WetArea(R,scaleF,Thres)
        WAR = computeWAR(R,Thres,scaleF);
        function [war] = computeWAR(m3d,thre,scaleF)
            war = nanmean(reshape((m3d>thre*scaleF),[],size(m3d,3)),1);
        end
    end
end
