filename = 'E:\OneDrive - Imperial College London\PHD IMPERIAL\PUBLICATIONS\2020-WaterLondon-\traced_maps\thames_water_wastewater_zones_traced.shp';
S = shaperead(filename);

for ID = 1:numel(S)
fill(S(ID).X([1:end-1,1]),S(ID).Y([1:end-1,1]),[0.5 0.5 0.5]);hold on;
end

UKMap = getUKMap();
plot(UKMap.borderE,UKMap.borderN,'k-');
axis equal

%%
filename = 'E:\OneDrive - Imperial College London\PHD IMPERIAL\PUBLICATIONS\2020-WaterLondon-\traced_maps\thames_water_wastewater_zones_traced.shp';
S = shaperead(filename);
region = struct;
region.minE = nanmin(arrayfun(@(s)nanmin(s.X)-1000,S));
region.maxE = nanmax(arrayfun(@(s)nanmax(s.X)+1000,S));
region.minN = nanmin(arrayfun(@(s)nanmin(s.Y)-1000,S));
region.maxN = nanmax(arrayfun(@(s)nanmax(s.Y)+1000,S));
unit = 1;
x_yr = region.minE*unit:1000:(region.maxE)*unit;
y_yr = region.minN*unit:1000:(region.maxN)*unit;
[XX,YY] = meshgrid(x_yr,y_yr);

%%
% XX = [523746];% Easting
% YY = [188202];% Northing

for YEAR = 2017:2019
    
    [DATA,status] = importNIMROD_P(XX,YY,YEAR);
    save(sprintf('PRS_London_%04d.mat',YEAR),'DATA','XX','YY','-v7.3');
    
end

time = datetime(2020,1,1):minutes(5):datetime(2020,10,4);
[DATA,status] = importNIMROD_P(XX,YY,time);
save(sprintf('PRS_London_%04d.mat',2020),'DATA','XX','YY','-v7.3');

%%
for YEAR = 2017:2020
    load(sprintf('PRS_London_%04d.mat',YEAR),'DATA','XX','YY');
    meanVal = nanmean(DATA,'spatial');
    pcolor(meanVal*24*365);shading flat
    cptcmap('precip_annualMeanUK', 'mapping','direct');%,'ncol',20);
    pause(0.5)
end

%%
filename = 'E:\OneDrive - Imperial College London\PHD IMPERIAL\PUBLICATIONS\2020-WaterLondon-\traced_maps\thames_water_wastewater_zones_traced.shp';
S = shaperead(filename);
region = struct;

DATA = [];
for YEAR = 2017:2020
    D = load(sprintf('PRS_London_%04d.mat',YEAR),'DATA','XX','YY');
    D.DATA = aggregate(D.DATA,1,hours(1)/minutes(5),'mm/h');
    DATA = appendTime(DATA,D.DATA);
end
rain = DATA.Val;
flag = rain>128*DATA.ScaleF;
Time = repmat(reshape(DATA.Time,[1,1,32929]),[size(DATA.Val,[1,2]),1]);
DATA.Val(flag) = DATA.Missingval;
%%
[subRain,dirRain,arealVal,flag] = deal([]);
for ID = 1:8
    [meanVal,subRain(:,ID),dirRain(:,ID),areaVal(:,ID),flag(:,ID)] = getSubArea(DATA,S(ID));
    pcolor(DATA.XX,DATA.YY,meanVal*24*365); shading flat
    cptcmap('precip_globalPrecip','mapping','scaled','ncol',1485,'flip',true);
    hold on
    plot(S(ID).X([1:end-1,1]),S(ID).Y([1:end-1,1]),'color','k');
    hold on;
end
totalRain = nansum(subRain.*areaVal,2)./nansum(areaVal);

ind = 1:32928;
totalRain = sum(subRain.*areaVal,2)./sum(areaVal);
totalDirRain = nanmean(buffer(totalRain,24),2);
DATA.Time.Format = 'default';
subarea = subRain(ind,:);
wholeLondon = totalRain(ind);
time = DATA.Time(ind);

subarea(isnan(subarea)) = -1;
wholeLondon(isnan(wholeLondon)) = -1;

A = table(subarea,wholeLondon,time);
writetable(A,'London2017_2020.csv');
A = readtable('London2017_2020.csv');
zone_name = [arrayfun(@(s)s.zone_name,S,'UniformOutput',false)]';
A.Properties.VariableNames(1:8) = zone_name;
writetable(A,'London2017_2020.csv');
%%

function [meanVal,subRain,dirRain,areaVal,flag] = getSubArea(DATA,shp)
isInArea = inpolygon(DATA.XX,DATA.YY,shp.X([1:end-1,1]),shp.Y([1:end-1,1]));

DATA.Val = reshape(DATA.Val,[],size(DATA.Val,3));
DATA.Val(~isInArea,:) = DATA.Missingval;
DATA.Val = reshape(DATA.Val,[size(DATA.XX),size(DATA.Val,2)]);
meanVal = nanmean(DATA,'spatial');
subRain = nanmean(DATA,'timeseries');
subRain = subRain';

dirRain = nanmean(buffer(subRain,24),2);

areaVal = nansum(~isnan(meanVal(:)));
flag = squeeze(nansum(DATA.Val>DATA.Missingval,[1,2]))./areaVal;

subRain(flag<=0.5)=NaN;
end


