% -----------------------------------------------------------------------%
% This script is to transform radar file from .gz to .mat for each day.
% Radar: NIMROD Composit Radar, from UK Met Office
% Resolution: 5min, 1km
% @ Yuting Chen
% Update: 2019.12
% -----------------------------------------------------------------------%

clear;clc
close all
cd('H:\CODE_MATLAB\SpatialTemporalDATA\Radar')
dataPath = 'K:\UkRadar2019\NIMROD_2016_2019\badc\ukmo-nimrod\data\composite\';
savePath = 'K:\UK_Radar_Matlab';


fid = fopen(['H:\CODE_MATLAB\SpatialTemporalDATA\RadarDialogFile\',...
    'NIMRODTRANSFORM_dialogFile.txt'],'w');


for year = 2019
    
    for mon = 1:12
        for day = 1%:eomday(year,mon)
            
            try
            tar_name = [dataPath,...
                sprintf(['uk-1km%s%04d%smetoffice-c-band-rain-radar_uk_%04d%02d%02d_1km-composite.dat.gz.tar'],...
                filesep,year,filesep,year,mon,day)];
            
            dirname = [savePath,'\backup'];
            delete([dirname,'\','*.gz']);% clean up folder befor untar.
            programFiles = untar(tar_name,dirname);
            
            pl = 0;
            tic
            sf = sprintf('%s%s%04d_%d_%d.mat',savePath,filesep,year,mon,day);
            [status, ~] = NIMROD_saveData(dirname,pl,sf);
            
            toc
            delete([dirname,'\','*.gz']);
            fprintf('Year%04dMon%02dDay%02d finished\n',year,mon,day);
            
            catch me
                fprintf(fid,'Check:%s\n',tar_name);
            end
            
        end
        fprintf('------Year%04dMon%02d finished-----\n',year,mon);
    end
end
