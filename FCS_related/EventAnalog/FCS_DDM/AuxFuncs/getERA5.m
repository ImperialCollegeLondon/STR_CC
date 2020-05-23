
function [CLIMATE,time,missval,dt] = getERA5(datainfo)

fileN = datainfo.fileN;
fileN_singleLevel = datainfo.fileN_singleLevel;
loni = datainfo.loni;
lati = datainfo.lati;
timei = datainfo.timestart;
lonlen = datainfo.lonlen;
latlen = datainfo.latlen;
timelen = datainfo.t_len;


getOneVar = @(fileN,var)ncread(fileN,var,[loni,lati,timei],[lonlen,latlen,timelen]);
GEOZ = getOneVar(fileN,'z');% 'geopotential';'m**2 s**-2';
O3 = getOneVar(fileN,'o3');%'Ozone mass mixing ratio' 'kg kg**-1'
RH = getOneVar(fileN,'r');% 'relative_humidity' '%';
CRWC = getOneVar(fileN,'crwc');% 'Specific rain water content' 'kg kg**-1'
T = getOneVar(fileN,'t');% 'air temperature';'k'
DT = getOneVar(fileN_singleLevel,'d2m');% '2 metre dewpoint temperature';'k'
MSP = getOneVar(fileN_singleLevel,'msl');%'Mean sea level pressure';'Pa'

% wind speed: 500hPa: for motion of many tropical cyclones
% wind speed: 700hPa: for shallower tropical cyclones
U = getOneVar(fileN,'u');% 'U component of wind' 'eastward_wind' 'm s**-1'
V = getOneVar(fileN,'v');% 'northward_wind' 'V component of wind' 'm s**-1'
time = datetime(datevec(ncread(fileN,'time',timei,timelen)/24+datenum(datetime(1900,1,1,0,0,0))));
missval = -32767;

CLIMATE = struct('GEOZ',GEOZ,'O3',O3,'RH',RH,'CRWC',CRWC,'T',T,'DT',DT,'U',U,'V',V,'MSP',MSP);

dt = 1;

end

