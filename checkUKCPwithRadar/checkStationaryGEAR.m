% ------------------------------------------------------------------------------- %
% This script is to Examine Stationary of time series of GEAR dataset.
% Time period: 1980-2000, 2007-2017
%
% @ Yuting Chen
% Imperial College London
% yuting.chen17@imperial.ac.uk
% Update: 2020.01.10
% ------------------------------------------------------------------------------- %

%% get Whole Time series for one location from GEAR
% REGION: whole UK
% resolution: 1 Km
% period: 1980-2017
% checking tips: ...


% Configuration
year0 = 1980;
year1 = 2017;

source = ['K:\GEAR\CEH_GEAR_daily_GB_',num2str(2010),'.nc'];
finfo = ncinfo(source);
varname = 'x';
x = ncread(source,varname); % easting-OSGB36 Grid reference
varname = 'y';
y = ncread(source,varname); % northing-OSGB36 Grid reference

RES = [];
tag = 1;
for bat = 387:416
    try
        Y_coor = y(3*bat+1 : 3*bat+3);
        X_coor = x(1:end);
        
        [X_coor,Y_coor] = meshgrid(X_coor,Y_coor);
        
        % ImportTimeseries
        RAIN = NaN(size(X_coor,1),size(X_coor,2),0);
        for year = year0:year1
            date0 = datenum(datetime(year,1,1):datetime(year,12,31));
            [XX,YY,R] = import_GEAR_DAILY(X_coor,Y_coor,date0,[]);
            RAIN = cat(3,RAIN,R);
        end
        % RAIN = squeeze(RAIN);
        
        % check Time series - Statistical Test
        for i = 1:size(RAIN,1)
            for j = 1:size(RAIN,2)
                [allTest{tag}] = checkST(squeeze(RAIN(i,j,:)));
                header = sprintf('----------Loc [%d, %d]---------',XX(i,j),YY(i,j));
                printTest(allTest{tag},header);
                tag = tag+1;
            end
        end
        fprintf('------------------------------------------\n');
        fprintf('Bat - %03d - finished. \n', bat);
        fprintf('------------------------------------------\n');
    catch
        1;
    end
end




% AUXILLARY FUNCTION
function [allTest] = checkST(RAIN)

% Detect Trend of TS

allTest = struct;

% Implement Baxter-King band-pass filter + Test
freq_min = 28;
freq_max = 31;
K = 12; % suggested by Baxter and King

try
    Y  = BK(RAIN,freq_min,freq_max,K);
    [allTest.adfH, allTest.adfpValue, stat, cValue, reg] = adftest(Y,'model','TS');
    
    Y = RAIN;
    [allTest.MKendallH, allTest.MKendallpValue] = Mann_Kendall_Modified(Y, 0.05);
    
    [allTest.SRhoTd, allTest.SRhopValue] = SpearmanRho(aggregate(Y,365,'mean'), 0.05);
catch
    [allTest.adfH, allTest.adfpValue, ...
        allTest.MKendallH, allTest.MKendallpValue, ...
        allTest.SRhoTd, allTest.SRhopValue]= deal(NaN);
end
end

function printTest(allTest,header)
if ~isnan(allTest.adfH)
    fprintf('%s\n',header);
    fprintf('----------------Test Result---------------\n')
    fprintf('adfTest\t %d \t pValue %f \n',allTest.adfH,allTest.adfpValue);
    fprintf('MKendall\t %d \t pValue %f \n',allTest.MKendallH,allTest.MKendallpValue);
    fprintf('SpearmanRho\t %d \t pValue %f \n',allTest.SRhoTd,allTest.SRhopValue);
    % fprintf('----------Test Result---------\n')
else
    1;
end
end



% % % load('H:\CODE_MATLAB\SpatialTemporalDATA\FloodMap\Birm.mat','Birm')
% % X_coor = 250269:1000:255000;%Loch Ness; %Birm.Easting;
% % Y_coor = 823116:1000:826000;%Loch Ness; %Birm.Northing;
% [X_coor,Y_coor] = meshgrid(X_coor,Y_coor);


