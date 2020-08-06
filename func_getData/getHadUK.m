function [HadUKMonth,E,N,HadUKAnnual] = getHadUK(var)
%
% @ Yuting
% Data source:
% http://data.ceda.ac.uk/badc/ukmo-hadobs/data/insitu/MOHC/HadOBS/HadUK-Grid/v1.0.1.0/1km/.../ann-20y/v20190808

arguments
    var (1,:) char = 'pr'
end

switch(var)
    case 'pr'
        varname = 'rainfall';
        filePath = 'K:\HadUK-Grid\';
        fileName_mon = 'rainfall_hadukgrid_uk_1km_mon-20y_198101-200012.nc';
        fileName_ann = 'rainfall_hadukgrid_uk_1km_ann-20y_198101-200012.nc';
    case 'tas'
        varname = 'tas';
        filePath = 'K:\HadUK-Grid\';
        fileName_mon = 'tas_hadukgrid_uk_1km_mon-20y_198101-200012.nc';
        fileName_ann = 'tas_hadukgrid_uk_1km_ann-20y_198101-200012.nc';
    otherwise
        error('Check input variable name.');
end

%% 20yr monthly
ncFileName = [filePath,fileName_mon];
S = ncinfo(ncFileName);

LAT=ncread(ncFileName,'latitude');
LON=ncread(ncFileName,'longitude');
[E, N] = ll2os(LAT, LON);
E = E/1000;
N = N/1000;

HadUKMonth = ncread(ncFileName,varname,[1,1,1],[Inf,Inf,Inf]);
HadUKMonth(HadUKMonth>1e10) = NaN;

%% 20yr
ncFileName = [filePath,fileName_ann];
HadUKAnnual = ncread(ncFileName,varname,[1,1,1],[Inf,Inf,Inf]);
HadUKAnnual(HadUKAnnual>1e10) = NaN;

end