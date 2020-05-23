function [HadUK,E,N,HadUKAnnual] = getHadUK()
%
% @ Yuting
% Data source: 
% http://data.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.0.1.0/1km/rainfall/ann-20y/v20190808

%% 20yr monthly
filePath = 'K:\HadUK-Grid\';
fileName = 'rainfall_hadukgrid_uk_1km_mon-20y_198101-200012.nc';
% 'rainfall_hadukgrid_uk_1km_ann-20y_198101-200012.nc';

ncFileName = [filePath,fileName];
S = ncinfo(ncFileName);

LAT=ncread(ncFileName,'latitude');
LON=ncread(ncFileName,'longitude');
[E, N] = ll2os(LAT, LON);
E = E/1000;
N = N/1000;

HadUK = ncread(ncFileName,'rainfall',...
    [1,1,1],...
    [Inf,Inf,Inf]);

HadUK(HadUK>1e10) = NaN;

%% 20yr
fileName = 'rainfall_hadukgrid_uk_1km_ann-20y_198101-200012.nc';
ncFileName = [filePath,fileName];
HadUKAnnual = ncread(ncFileName,'rainfall',...
    [1,1,1],...
    [Inf,Inf,Inf]);
HadUKAnnual(HadUKAnnual>1e10) = NaN;

end