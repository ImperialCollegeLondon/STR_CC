
clear;clc

% update chunksize for more effecient .nc
try
    
    for year = 2011:2018
        for mon = 1:12
            
            tic
            
            TS = datetime(year,mon,1,0,0,0):1/24:datetime(year,mon,eomday(year,mon),23,0,0);
            time = datenum(TS);
            
            fileName = sprintf('pr_nimrod_uk_2.2km_1hr_%04d%02d.nc',year,mon);
            
            filePath = 'K:\UK_Radar_NetCDF\';
            
            pr = ncread([filePath,fileName],'pr');
            E = ncread([filePath,fileName],'E');
            N = ncread([filePath,fileName],'N');
            
            sfPath = 'K:\UK_Radar_NetCDF_good\';
            createRadarNetCDF(sfPath,fileName,N,E,time,pr);
            
            toc
            
        end
    end
    
catch me
    
    1;
    
end
