function Config = getConfig(regionname,MON,Period,ENSEMBLENO)

if iscell(ENSEMBLENO) && length(ENSEMBLENO) == 1
    ENSEMBLENO = ENSEMBLENO{1};
end
REGIONS = REGIONS_info();
region = getfield(REGIONS,regionname); %#ok<GFLD>
IntFac = 32; % to make the resolution = (Nimrod)
if strcmp(Period,'1980-2000')
    file = ['D:\UKCP18\',region.Name,'\',sprintf('Ensems_%s_mon%02d.mat',ENSEMBLENO,MON)];
elseif strcmp(Period,'2020-2040')
    file = ['D:\UKCP18_Future\',region.Name,'\',sprintf('Ensems_%s_mon%02d.mat',ENSEMBLENO,MON)];
elseif strcmp(Period,'2060-2080')
    file = ['D:\UKCP18_Future_2060_2080\',region.Name,'\',sprintf('Ensems_%s_mon%02d.mat',ENSEMBLENO,MON)];
elseif strcmp(Period,'2007-2018') % && strcmp(ENSEMBLENO,'RAD')
    file = ['D:\UKCP18\',region.Name,'\',sprintf('Radar_mon%02d.mat',MON)];
end

method = struct;
method.eventSeperation = 'WAR-based';
method.MCScrit = 'PIMF>5';
method.DURcrit = 'WAR>0.05';

data = struct;
data.file = file;
data.ENSNO = ENSEMBLENO;
data.IntFac = IntFac;

saveIt = struct;
saveIt.path = 'D:\UKCP18\ConvectiveStorms';

Config = struct('data',data,...
    'region',region,...
    'Month',MON,...
    'method',method,...
    'saveIt',saveIt);

end
