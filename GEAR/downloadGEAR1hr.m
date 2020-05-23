
YEARRANGE = [1990:2014];
options = weboptions('Username',getYutingEmail(),'Password',getYutingEmail('PAssword'),'Timeout',Inf);
weblink = 'https://catalogue.ceh.ac.uk/datastore/eidchub/d4ddc781-25f3-423a-bba0-747cc82dc6fa//';


for year = YEARRANGE
    for mon = 1:12
        fullFilename = [weblink,sprintf('CEH-GEAR-1hr_%04d%02d.nc',year,mon)];
        OFN = websave(['K:\GEAR-1hr\',sprintf('CEH-GEAR-1hr_%04d%02d.nc',year,mon)],...
            fullFilename,options);
        fprintf('Finished:CEH-GEAR-1hr_%04d%02d.nc\n',year,mon);
    end
    
end

%CEH-GEAR-1hr_199001.nc