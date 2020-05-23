% extract radar data

XX = [523746];% Easting
YY = [188202];% Northing
% PRS = [];
% PTime = [];

for YEAR = 2019 % 2006:2019
    PRS = [];PTime = [];
    [DATA,status] = importNIMROD_P(XX,YY,YEAR);
    PTime = [PTime;DATA.Time];
    PRS = cat(1,PRS,squeeze(DATA.Val));
    save(sprintf('PRS_%04d.mat',YEAR),'PRS','PTime','-v7.3');
end

% %%
% Time = PTime;
% Time.Format = 'dd-MMM-uuuu HH:mm:SS';
% 
% Intensity = double(PRS)/32;
% Intensity(Intensity<0) = -999;
% %%
% Tab = table(Time,Intensity);
% writetable(Tab,'H:\DATA_CEDA\STATION_hourly_2006_2020\radar2018.csv');