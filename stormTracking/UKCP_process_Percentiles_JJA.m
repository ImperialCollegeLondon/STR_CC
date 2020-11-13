
clear;clc
close all
prcval = [97.5,99,99.9,99.99];
warning on

for regionName = {'CPM_NW'};%{'CPM_S','CPM_NE'}%{'EUK','SCO','WAL'},% % 'CPM_NE'
    regionName = regionName{1};
    PRC1 = getPRC('1980-2000',regionName);
    PRC2 = getPRC('2060-2080',regionName);
    save([regionName,'_prctile.mat'],'PRC1','PRC2','prcval');
end

function PRC = getPRC(Period,regionName)
PRC = [];
ENSEMBLENO = getEnsNos();
for ensNo = ENSEMBLENO([1:12])%{'RAD'}%
    ensNo = ensNo{1};%
    MON = [6:8];
    [prc,Config] = getData(regionName,MON,Period,ensNo);
    PRC = [PRC;prc]; 
end
fprintf('%s %s done\n',regionName,Period);
end

function [prc,Config] = getData(regionName,MON,Period,ENSNO)

RainEnsembles = [];

UKMap = getUKMap();
[E,N] = getEN(getfield(REGIONS_info(),regionName));
in = inpolygon(E,N,UKMap.borderE/1000,UKMap.borderN/1000);


for mon = MON(:)'
    Config = getConfig(upper(regionName),mon,Period,ENSNO);
    datafile = Config.data.file;
    A = load(datafile);RainEnsembles_temp = [];
    if isfield(A,'RainEnsembles')
        RainEnsembles_temp = A.RainEnsembles;
    else
        RainEnsembles_temp{1} = A.Rain;
    end
    clear A
    for ensNo = 1:length(RainEnsembles_temp)
        % # only valid for UKCP hourly projections #
        if iscell(RainEnsembles) && ~isempty(RainEnsembles{ensNo})
            RainEnsembles{ensNo} = cat(3,RainEnsembles{ensNo},...
                RainEnsembles_temp{ensNo});
        else
            RainEnsembles{ensNo} = RainEnsembles_temp{ensNo};
        end
    end
    % fprintf('Output is int R with scaleF\n');
end

fprintf('Output is int R with scaleF\n');
Config.Month = MON;

tic
In = in(:);
R = reshape(RainEnsembles{1},[],size(RainEnsembles{1},3));
R(~In,:) = NaN;
toc
R = reshape(R,[],1);
prc = prctile(R,[97.5,99,99.9,99.99]);
prc = double(prc)/Config.data.IntFac;
end
