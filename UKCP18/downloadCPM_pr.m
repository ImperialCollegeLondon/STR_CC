% ----------------------------------------------------------------------- %
% 
% DOWNLOADING FILE
% suitable for UKCP CPM data
% 
% ----------------------------------------------------------------------- %

clear all
close all
clc


ftp_client = ftp('ftp.ceda.ac.uk','ychen021','AaBb14207');


%% Download data

ENSEMBLENO = {'01','04','05','06','07','08','09','10','11','12','13','15'};


for M=1:length(ENSEMBLENO)
    
    filePath = ['J:/UKCP18_CPM_Future/',ENSEMBLENO{M},'/'];
    cd(ftp_client, sprintf('/badc/ukcp18/data/land-cpm/uk/2.2km/rcp85/%s/pr/1hr/latest/',ENSEMBLENO{M}));
    fol = sprintf('/badc/ukcp18/data/land-cpm/uk/2.2km/rcp85/%s/pr/1hr/latest/',ENSEMBLENO{M});
    
    for mon = [1:12]
        
        details = dir(ftp_client,fol);
        
        for year = [2060:2080] % 2020:2040% 1980:2000
            
            if ~((year==2060 && mon ~=12) || (year == 2080 && mon == 12)...
                    || (year==2020 && mon ~=12) || (year == 2040 && mon == 12))
                
                fileName = sprintf('pr_rcp85_land-cpm_uk_2.2km_%s_1hr_%04d%02d01-%04d%02d30.nc',...
                    ENSEMBLENO{M},year,mon,year,mon);
                listaRain = dir([filePath,fileName]);
                
                if isempty(listaRain)
                    try
                        mget(ftp_client,fileName,filePath);
                    catch me
                        try
                            mget(ftp_client,[fileName,'.gz'],filePath);
                            gunzip([filePath,'/',fileName,'.gz'], filePath);
                            delete([filePath,'/',fileName,'.gz']);
                        catch me
                            error('Check')
                        end
                    end
                else
                    % error('Check')
                end
            end
            fprintf('%s-Mon%02d-Year%04d finished scanning.\n',ENSEMBLENO{M},mon,year);
            clear mex
        end
    end
end

close (ftp_client)



