% Process whole UK Pr


%% CONFIGURATION
clear;clc

addpath(genpath('C:\Users\Yuting Chen\Dropbox (Personal)\MatLab_Func'));
addpath(genpath(cd));


REGIONS = REGIONS_info();
dataSP = 'H:\DATA_CLIMATE\UKCP18\';
region = REGIONS.REGIONS.London;%SWestuk;%Westuk;%Scotland;
Period = '1980-2000';

itag = 0;

LetItRun(dataSP,REGIONS.SWestuk,Period);

LetItRun(dataSP,REGIONS.Westuk,Period);

LetItRun(dataSP,REGIONS.Scotland,Period);

LetItRun(dataSP,REGIONS.London,Period);

fprintf('Finished\n')

function LetItRun(dataSP,region,Period)

[E,N] = getEN(region);

% Compute the GridAreaMatrix
region.gridArea = diff([E(1,1)-1.1;E(:,1)+1.1])...
    *diff([N(1,:)-1.1,N(1,end)+1.1]);

%% Part *: Identify Cells

% UKCP
THRE = [1 3 5 7 12 20 40 50 100];

for season = 1:4
    CELLA = cell(1,length(THRE));
    [E,N,RainEnsembles,scaleF,region] = loadRainEnsembles(region,getMons(season),Period);
    for threi = 1:length(THRE)
        thre = THRE(threi);
        tic
        for enNo = 1:length(RainEnsembles)
            rr = RainEnsembles{enNo};
            rr(rr<=thre*scaleF) = 0;
            rr = logical(rr);
            CELLA{threi}{enNo} = [];
            for ti = 1:size(rr,3)
                [output] = getCells(squeeze(rr(:,:,ti)));
                CELLA{threi}{enNo} = [CELLA{threi}{enNo};output];
            end
            fprintf('EnNo%02d,Thre%02d,Season%02d Done----\n',enNo,thre,season);
        end
        
        toc
        fprintf('----Thre%02d,Season%02d Done----\n',thre,season);
    end
    save(sprintf('%sCellProp_UKCP_Season%01d_%s.mat',dataSP,season,region.Name),'CELLA','THRE');
end
%%

% NIMROD
for season = 1:4
    CELLA = cell(1,length(THRE));
    [E,N,RainEnsembles,scaleF,region] = loadRainNimrod(region,getMons(season));
    RainEnsembles = {RainEnsembles};
    for threi = 1:length(THRE)
        thre = THRE(threi);
        tic
        for enNo = 1:length(RainEnsembles)
            rr = RainEnsembles{enNo};
            rr(rr<=thre*scaleF) = 0;
            rr = logical(rr);
            CELLA{threi}{enNo} = [];
            for ti = 1:size(rr,3)
                [output] = getCells(squeeze(rr(:,:,ti)));
                CELLA{threi}{enNo} = [CELLA{threi}{enNo};output];
            end
            fprintf('EnNo%02d,Thre%02d,Season%02d Done----\n',enNo,thre,season);
        end
        
        toc
        fprintf('----Thre%02d,Season%02d Done----\n',thre,season);
    end
    save(sprintf('%sCellProp_NIMROD_Season%01d_%s.mat',dataSP,season,region.Name),...
        'CELLA','THRE');
end


    function [output] = getCells(BW)
        
        % BW = logical(rr);
        % [B,L] = bwboundaries(BW,'noholes');
        stats = regionprops('Table',BW,'Area');
        % 'Centroid','MajorAxisLength','MinorAxisLength',
        % imshow(label2rgb(L, @jet, [.5 .5 .5]))
        output = stats.Area;
        
    end


end
