% Example:
% time = 1:720;
% N = 1:606;
% E = 1:484;
% pr = rand(length(E),length(N),length(time));
% 
% filePath = [cd,'\'];
% fileName = 'XXXX.nc';
% createRadarNetCDF(filePath,fileName,N,E,time,pr)
%
% scaleF: 1
% data format: single

function createRadarNetCDF(filePath,fileName,N,E,time,pr)
%Open the file
ncid = netcdf.create([filePath,fileName],'NETCDF4');

%Define the dimensions
dimidE = netcdf.defDim(ncid,'E',length(E));
dimidN = netcdf.defDim(ncid,'N',length(N));
dimidt = netcdf.defDim(ncid,'time',length(time));

%Define IDs for the dimension variables (pressure,time,latitude,...)
E_ID=netcdf.defVar(ncid,'E','NC_DOUBLE',[dimidE]);
N_ID=netcdf.defVar(ncid,'N','NC_DOUBLE',[dimidN]);
date_ID=netcdf.defVar(ncid,'time','NC_DOUBLE',[dimidt]);

%Define the main variable (pr)
pr_ID = netcdf.defVar(ncid,'pr','NC_FLOAT',[dimidE dimidN dimidt]);% single format as MetUM.
netcdf.defVarChunking(ncid,pr_ID,'CHUNKED',[length(E) length(N) 1]);

%We are done defining the NetCdf
netcdf.endDef(ncid);

%Then store the dimension variables in
netcdf.putVar(ncid,E_ID,E);
netcdf.putVar(ncid,N_ID,N);
netcdf.putVar(ncid,date_ID,time);

%Then store my main variable
netcdf.putVar(ncid,pr_ID,pr);

%We're done, close the netcdf
netcdf.close(ncid);

clear mex

end

