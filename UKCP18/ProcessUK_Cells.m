% Process regional Pr
% London Region;

%% CONFIGURATION
clear;clc

addpath(genpath('C:\Users\Yuting Chen\Dropbox (Personal)\MatLab_Func'));
addpath(genpath(cd));


REGIONS = REGIONS_info();
dataSP = 'H:\DATA_CLIMATE\UKCP18\';


LetItRun(dataSP,REGIONS.UK);


fprintf('Finished\n')

function LetItRun(dataSP,region)


%% Part *: Identify Cells
%{

% UKCP
THRE = [1,3,20];%[1 3 5 7 12 20 40 50 100];
ENSEMBLENO=getEnsNos();

for season = 1:4
    for enNo = 1:length(ENSEMBLENO)
        for mon = getMons(season)
            tic
            CELLA = cell(1,length(THRE));
            for year  = 1980:2000%% to be customized due to histoData staring from 198012
                if ~((year == 1980 && mon ~= 12) || (year == 2000 && mon == 12))
                    
                    
                    tag = 0;
                    imageNo = 1;
                    time0 = datetime(year,mon,1,0,0,0);
                    while (tag == 0)
                        try
                            [E,N,RainEnsembles,scaleF,region] = readCPM_nc_1(region,year,mon,ENSEMBLENO{enNo},imageNo);
                        catch me
                            tag = 1;
                        end
                        
                        timeNum = datenum(time0+hours(imageNo-1));
                        for threi = 1:length(THRE)
                            rr = RainEnsembles;
                            thre = THRE(threi);
                            rr(rr<=thre*scaleF) = 0;
                            rr = logical(rr);
                            [output] = getCells(rr,timeNum);
                            CELLA{threi}= [CELLA{threi};output];
                        end
                        imageNo = imageNo+1;
                    end
                    fprintf('----%s-%04dSeason%02dMon%02d Done----\n',...
                        ENSEMBLENO{enNo},year,season,mon);
                    
                end
            end
            
            save(sprintf('%sCellProp_UKCP%s_Mon%01d_%s.mat',...
                dataSP,ENSEMBLENO{enNo},mon,region.Name),'CELLA','THRE');
            toc
            
        end
    end
end

%}


%%
% NIMROD
THRE = [1,3,20];%[1 3 5 7 12 20 40 50 100];

for season = 1:4
        for mon = getMons(season)
            tic
            CELLA = cell(1,length(THRE));
            for year  = 2007:2018%% to be customized due to histoData staring from 198012

                    tag = 0;
                    imageNo = 1;
                    time0 = datetime(year,mon,1,0,0,0);
                    while (tag == 0)
                        try
                            [E,N,RainObs,scaleF,region] = readRAD_nc_1(region,year,mon,imageNo);
                        catch me
                            tag = 1;
                        end
                        rr = RainObs;
                        timeNum = datenum(time0+hours(imageNo-1));
                        for threi = 1:length(THRE)
                            rr = RainObs;
                            thre = THRE(threi);
                            rr(rr<=thre*scaleF | isnan(rr)) = 0;
                            rr = logical(rr);
                            [output] = getCells(rr,timeNum);
                            CELLA{threi}= [CELLA{threi};output];
                        end
                        imageNo = imageNo+1;
                    end
                    fprintf('----%04dSeason%02dMon%02d Done----\n',...
                        year,season,mon);

            end
            
            save(sprintf('%sCellProp_RAD_Mon%01d_%s.mat',...
                dataSP,mon,region.Name),'CELLA','THRE');
            toc
            
        end

end


    function [output] = getCells(BW,timeNum)
        
        % BW = logical(rr);
        % [B,L] = bwboundaries(BW,'noholes');
        stats = regionprops('Table',BW,'Area','Centroid','Eccentricity');
        stats.time = timeNum.*ones(size(stats,1),1);
        % 'Centroid','MajorAxisLength','MinorAxisLength',
        % imshow(label2rgb(L, @jet, [.5 .5 .5]))
        output = stats;%stats.Area;
        
    end


end

%
% read RAD one snapshot.
function [E,N,Rain,scaleF,region] = readRAD_nc_1(region,year,mon,imageNo);

fileName = sprintf('pr_nimrod_uk_2.2km_1hr_%04d%02d.nc',...
    year,mon);
filePath = ['K:/UK_Radar_NetCDF/'];
listaRain = [filePath,fileName];
A = ncinfo(listaRain);
% LAT=ncread(listaRain,'latitude');
% LON=ncread(listaRain,'longitude');
Rain = squeeze(ncread(listaRain,'pr',...
    [1,    1,  imageNo],...
    [Inf,Inf,    1]));
scaleF = 1;

E = NaN;
N = NaN;% can be acquired using function in 'getAllCellsUK.m'
end



% read CPM one snapshot.
function [E,N,Rain,scaleF,region] = readCPM_nc_1(region,year,mon,ensNo,imageNo);

fileName = sprintf('pr_rcp85_land-cpm_uk_2.2km_%s_1hr_%04d%02d01-%04d%02d30.nc',...
    ensNo,year,mon,year,mon);
filePath = ['K:/UkCp18/',ensNo,'/'];
listaRain = [filePath,fileName];
A = ncinfo(listaRain);
% LAT=ncread(listaRain,'latitude');
% LON=ncread(listaRain,'longitude');
Rain = squeeze(ncread(listaRain,'pr',...
    [1,    1,  imageNo,  1],...
    [Inf,Inf,    1,  1]));
scaleF = 1;

E = NaN;
N = NaN;% can be acquired using function in 'getAllCellsUK.m'
end
