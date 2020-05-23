%--------------------------------------------------%
% this file is to read monthly meanensembles of CPM Pr
%
% @ Yuting Chen
% Imperial College London
% yuting.chen17@imperial.ac.uk
% Update: 2019.12.**
%--------------------------------------------------%

function [RainEnsembles,IntFac,E,N] = readCPM_pr_average(region,ENSEMBLENO,MON,IntFac)
%
% READCPM_PR() gives .....
%
% Example:
%         IntFac = 1;
%         options = 'Mean';
%         [RainEnsembles,IntFac,E,N] = readCPM_pr(region,ENSEMBLENO,MON,IntFac);
%
%
%
% @ Yuting Chen
% Imperial College London


arguments
    
    region (1,1) struct
    ENSEMBLENO (1,:) cell
    MON (1,1) double
    IntFac (1,1) double = 1;
    
end


%% Configuration
filePath = ['K:/UkCp18/',ENSEMBLENO{1},'/pr_rcp85_land-cpm_uk_2.2km_'];
ncFileName = [filePath,ENSEMBLENO{1},'_1hr_19910801-19910830.nc'];
LAT=ncread(ncFileName,'latitude');
LON=ncread(ncFileName,'longitude');
[E, N] = ll2os(LAT, LON);
E = E/1000;
N = N/1000;

%%
[region.i,region.j] = getRegionIJ(E,N,region.minE,region.minN);

warning off

for mon = MON
    
    RainEnsembles = cell(length(ENSEMBLENO),1);
    
    for M=1:length(ENSEMBLENO)
        
        filePath = ['K:/UkCp18/',ENSEMBLENO{M},'/'];
        
        MRain = [];
        
        L = 1;
        
        for year = 1980:2000
            tic
            if ~((year==1980 && mon ~=12) || (year == 2000 && mon == 12))
                
                fileName = sprintf('pr_rcp85_land-cpm_uk_2.2km_%s_1hr_%04d%02d01-%04d%02d30.nc',...
                    ENSEMBLENO{M},year,mon,year,mon);
                
                listaRain = dir([filePath,fileName]);
                
                
                if ~isempty(listaRain)
                    
                    listaRain = fullfile({listaRain.folder},{listaRain.name});
                    listaRain = listaRain{1};
                    A = ncinfo(listaRain);
                    
                    MRain(:,:,L) = nanmean(squeeze(ncread(listaRain,'pr',...
                        [1,1,1,1],...
                        [A.Variables(1).Dimensions(:,1).Length,...
                        A.Variables(1).Dimensions(:,2).Length,...
                        A.Variables(1).Dimensions(:,3).Length,...
                        A.Variables(1).Dimensions(:,4).Length])),3);
                    
                    
                else
                    
                end
                L = L+1;
                
            end
            toc
        end
        
        RainEnsembles{M} = MRain;
        fprintf('Ensembles %02d\n',M);
        
    end
    
    
    IntFac = 1;
    fprintf('Mon %02d Exported out, in <double> format\n',M);
    
    
    fprintf('---%s-Mon%02d---\n',region.Name,mon);
    
end



filePath = ['K:/UkCp18/01/pr_rcp85_land-cpm_uk_2.2km_'];
ncFileName = [filePath,'01','_1hr_19910801-19910830.nc'];
LAT=ncread(ncFileName,'latitude');
LON=ncread(ncFileName,'longitude');

[E, N] = ll2os(LAT, LON);
E = E/1000;
N = N/1000;


end

