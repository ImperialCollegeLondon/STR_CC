clear;
clc;
addpath(genpath(cd))
addpath(genpath('C:\Users\Yuting Chen\Dropbox (Personal)\MatLab_Func'));

%% Configuration
ENSEMBLENO=getEnsNos();
filePath = ['J:/UkCp18/',ENSEMBLENO{1},'/pr_rcp85_land-cpm_uk_2.2km_'];
ncFileName = [filePath,ENSEMBLENO{1},'_1hr_19910801-19910830.nc'];
LAT=ncread(ncFileName,'latitude');
LON=ncread(ncFileName,'longitude');
[E, N] = ll2os(LAT, LON);
E = E/1000;
N = N/1000;

UKMap = load('H:\CODE_MATLAB\SpatialTemporalDATA\shapeFileFolder\UKBorderGrid.mat');
warning off
MERain = [];%[mon,enNo,[2d matrix]]

for mon = 1:12
    RainEnsembles = cell(length(ENSEMBLENO),1);

    for M=1:length(ENSEMBLENO)
        
        filePath = ['J:/UkCp18/',ENSEMBLENO{M},'/pr_rcp85_land-cpm_uk_2.2km_'];
        
        listaRain=dir([filePath,ENSEMBLENO{M},'_1hr_*',sprintf('%02d',mon),'30.nc']);
        
        listaRain=fullfile({listaRain.folder},{listaRain.name});

        RR = [];
        for L=1:length(listaRain)
            
            A=ncinfo(listaRain{L});
            
            rr = ncread(listaRain{L},'pr',...
                [1,1,1,1],...
                [Inf,Inf,Inf,...
                A.Variables(1).Dimensions(:,4).Length]);
            % rr(rr<=1/32) = NaN;
            RR(L,:,:) = squeeze(nanmean(rr,3));
        end
        rr = squeeze(nanmean(RR,1));
        MERain(mon,M,:,:) = rr;%[mon,enNo,[2d matrix]]
        fprintf('---Mon%02d Ensembles%02d---\n',mon,M);
    end
    fprintf('---%s-Mon%02d---\n','Whole UK',mon);
end
save('UKCPRainPattern.mat','MERain');




