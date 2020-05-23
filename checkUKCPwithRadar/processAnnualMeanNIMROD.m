% ------------------------------------------------------------------------%
% This file will:
% Compute Annual mean from NIMROD Radar.
% It will be done by computing the average rainfall pattern for each month.
% Available period of NIMROD Radar: 2007-2018.
% Information will be saved in 'D:/UKCP18/NIMRODRainPattern.mat';
%
% @ Yuting Chen
% Update: 2019.12.11
% ------------------------------------------------------------------------%

clear;
clc;
addpath(genpath(cd))
addpath(genpath('C:\Users\Yuting Chen\Dropbox (Personal)\MatLab_Func'));

YEAR = 2007:2018;
REGIONS = REGIONS_info();
region = REGIONS.UK;

%% Configuration

MERain = [];%[mon,yearNo,[2d matrix]]

for mon = 1:12
    
    for yearNo = 1:length(YEAR)
        
        [E,N,Rain,scaleF,region] = readRAD_nc_month(region,YEAR(yearNo),mon);
        
        rr = squeeze(nanmean(double(Rain)/scaleF,3));
        
        MERain(mon,yearNo,:,:) = rr;
        
        fprintf('---Mon%02d Year%4d---\n',mon,YEAR(yearNo));
        
    end
    
    fprintf('---%s-Celendar Month %02d Finished---\n','Whole UK',mon);
    
end

save('D:/UKCP18/NIMRODRainPattern.mat','MERain','E','N','YEAR','-v7.3');



%% AUXILLAY FUNCITON
function [E,N,Rain,scaleF,region] = readRAD_nc_month(region,year,mon);
%
% read RAD.nc for $year$.$mon$.$region$
% resolution: 1hour;
%

fileName = sprintf('pr_nimrod_uk_2.2km_1hr_%04d%02d.nc',...
    year,mon);
filePath = ['K:/UK_Radar_NetCDF/'];
listaRain = [filePath,fileName];

A = ncinfo(listaRain);

% LAT = ncread(listaRain,'latitude');
% LON = ncread(listaRain,'longitude');

tic
Rain = squeeze(ncread(listaRain,'pr',...
    [1,    1,  1],...
    [region.dimE,region.dimN,eomday(year,mon)*24]));
scaleF = 1;
toc

clear mex

E = region.minE:2.2:region.minE+(region.dimE-1)*2.2;
N = region.minN:2.2:region.minN+(region.dimN-1)*2.2;
[N,E] = meshgrid(N,E);
end





