% extract radar data
clear;clc

REGIONS = REGIONS_info();
region = REGIONS.Westuk;

%% BIRM 110Km*110Km
unit = 1000;
x_yr = region.minE*unit:1000:(region.minE+109)*unit;
y_yr = region.minN*unit:1000:(region.minN+109)*unit;
[XX,YY] = meshgrid(x_yr,y_yr);
%%
% XX = [523746];% Easting
% YY = [188202];% Northing

for YEAR = 2013 %[2006,2012,2013]
    
    [DATA,status] = importNIMROD_P(XX,YY,YEAR);
    save(sprintf('PRS_Birm_%04d.mat',YEAR),'DATA','-v7.3');
    
end

%% BIRM 500Km*500Km
unit = 1000;
XDim = 250;%500;
x_yr = (region.minE-(XDim/2-55))*unit:1000:(region.minE+(XDim/2+54))*unit;
y_yr = (region.minN-(XDim/2-55))*unit:1000:(region.minN+(XDim/2+54))*unit;
[XX,YY] = meshgrid(x_yr,y_yr);

%%
% XX = [523746];% Easting
% YY = [188202];% Northing

for YEAR = 2013 %[2006,2012,2013]
    PRS = [];PTime = [];
    [DATA,status] = importNIMROD_P(XX,YY,YEAR);
    
    save(sprintf('H:\CODE_MATLAB\PRS_Birm250_%04d.mat',YEAR),'DATA','-v7.3');
end




% %%
% Time = PTime';
% Time.Format = 'dd-MMM-uuuu HH:mm:SS';
% 
% Intensity = double(PRS)/32;
% Intensity(Intensity<0) = -999;
% %%
% Tab = table(Time,Intensity);
% writetable(Tab,'H:\DATA_CEDA\STATION_hourly_2006_2020\radar2018.csv');