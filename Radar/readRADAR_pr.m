%--------------------------------------------------%
%
% this file is to read ensembles of CPM Pr.
% file is saved in .mat for quicker load next time.
% a factor of 32 was used to squeeze the filesize
% Data output:
%    Rain <int16> [E,N,T]
%    both E and N dimension: are in ascending order
% 
% @ Yuting Chen
% Imperial College London
% yuting.chen17@imperial.ac.uk
% Update: 2019.12.**
%--------------------------------------------------%


clear;
clc;


%% Configuration

IntFac = 32; % to make the resolution = (Nimrod)

REGIONS = REGIONS_info();

letItRun(REGIONS.SCO,IntFac)
letItRun(REGIONS.EUK,IntFac)
letItRun(REGIONS.WAL,IntFac)

letItRun(REGIONS.SWestuk,IntFac)
letItRun(REGIONS.Westuk,IntFac)
letItRun(REGIONS.Scotland,IntFac)
letItRun(REGIONS.London,IntFac)

function letItRun(region,IntFac)

%% Let it Run 
filePath = ['K:/UK_Radar_NetCDF/'];
ncFileName = [filePath,'pr_nimrod_uk_2.2km_1hr_201809.nc'];
E = ncread(ncFileName,'E');
N = ncread(ncFileName,'N');


%%
% to do:
% delete ensembleNo related

%%

[region.i,region.j] = getRegionIJ(E,N,region.minE,region.minN);

warning off

L = 1;
for mon = 6:8%1:12
    
    
    
    mkdir(['D:/UKCP18/',region.Name]);
    saveDir = ['D:/UKCP18/',region.Name,'/Radar',...
        '_mon',sprintf('%02d',mon)];
    
    Rain = int16(NaN(region.dimE,region.dimN,11*30*24));%preallocate space
    
    L = 1;
    
    for year = 2007:2018
            
            fileName = sprintf('pr_nimrod_uk_2.2km_1hr_%04d%02d.nc',...
                year,mon);
            
            listaRain = [filePath,fileName];
            
            A = ncinfo(listaRain);
            
            Rain(:,:,L:L+eomday(year,mon)*24-1) = round(IntFac*squeeze(ncread(listaRain,'pr',...
                [region.i,region.j,1],...
                [region.dimE,region.dimN,A.Dimensions(3).Length])));

            L = L+eomday(year,mon)*24;
    end
    

    if getfield(whos('Rain'),'bytes') > 5e8
        uiwait(msgbox('Almost out of memory, data will be saved in a new file',...
            'Oooops','__'));
    end
    
    
    save([saveDir,'.mat'],'Rain','-v7.3');
    
    fprintf('---%s-Mon%02d---\n',region.Name,mon);
end

end



% load('D:\UKCP18\SWestuk\Radar_mon01.mat');
% for i = 1:size(Rain,3)
%     
%     pcolor(squeeze(double(Rain(:,:,i)))/32);
%     shading flat
%     caxis([0,10])
%     cptcmap('GMT_drywet', 'mapping', 'scaled');
%     colorbar
%     
%     pause(0.1);
% end
