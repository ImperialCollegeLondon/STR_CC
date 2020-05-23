clear;clc
clear mex

filePath = 'D:/UKCP18/';
fileName = 'NIMRODRainPattern';
load([filePath,fileName,'.mat'],'MERain');%[Mon,Year,2D Matrix]

REGIONS = REGIONS_info();
region = REGIONS.UK;

E = region.minE:2.2:region.minE+(region.dimE-1)*2.2;
N = region.minN:2.2:region.minN+(region.dimN-1)*2.2;
% [N,E] = meshgrid(N,E);% not same as X-coor,Y-coor

try
    pr = permute(MERain,[3,4,1,2]);
    createRadarNetCDF_MonYearLoc(filePath,[fileName,'.nc'],N,E,1:12,2007:2018,pr)
catch me
end

function createRadarNetCDF_MonYearLoc(filePath,fileName,N,E,Month,Year,pr)
%Open the file
ncid = netcdf.create([filePath,fileName],'NETCDF4');

%Define the dimensions
dimidE = netcdf.defDim(ncid,'E',length(E));
dimidN = netcdf.defDim(ncid,'N',length(N));
dimidt = netcdf.defDim(ncid,'Month',length(Month));
dimidy = netcdf.defDim(ncid,'Year',length(Year));

%Define IDs for the dimension variables (pressure,time,latitude,...)
E_ID = netcdf.defVar(ncid,'E','NC_DOUBLE',[dimidE]);
N_ID = netcdf.defVar(ncid,'N','NC_DOUBLE',[dimidN]);
month_ID = netcdf.defVar(ncid,'Month','NC_DOUBLE',[dimidt]);
year_ID = netcdf.defVar(ncid,'Year','NC_DOUBLE',[dimidy]);

%Define the main variable (pr)
pr_ID = netcdf.defVar(ncid,'pr','NC_FLOAT',[dimidE dimidN dimidt dimidy]);% <single>, same as MetUM.
netcdf.defVarChunking(ncid,pr_ID,'CHUNKED',[length(E) length(N) 1 1]);

%We are done defining the NetCdf
netcdf.endDef(ncid);

%Then store the dimension variables in
netcdf.putVar(ncid,E_ID,E);
netcdf.putVar(ncid,N_ID,N);
netcdf.putVar(ncid,month_ID,Month);
netcdf.putVar(ncid,year_ID,Year);

%Then store my main variable
netcdf.putVar(ncid,pr_ID,pr);

%We're done, close the netcdf
netcdf.close(ncid);

clear mex

end

















