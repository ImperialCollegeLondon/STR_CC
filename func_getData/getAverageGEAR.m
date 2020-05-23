
function [GEAR,X_coor,Y_coor,DIST] = getAverageGEAR(region,scale,YEARVEC)
% GETAVERAGEGEAR(region,scale) give mean pr for each year, default range is 
% <1980:1999>
% area: done for whole UK
%
% Output format: 
%               GEAR: 3D # <double> [yearNo/monthNo,N,E]
%               X_coor: 2D # unit:km
%               Y_coor: 2D # unit:km
% Example:
%         [GEAR,X_coor,Y_coor] = getAverageGEAR(~,'year')
%         or
%         [GEAR,X_coor,Y_coor] = getAverageGEAR(~,'month')
%         or
%         [GEAR,X_coor,Y_coor] = getAverageGEAR(~,'year',[1990:1995])
%
% @ Yt
%
arguments
    
    region (1,1) struct
    scale (1,:) string
    YEARVEC (1,:) double = 1981:2000

end


source = ['K:\GEAR-month\CEH_GEAR_monthly_GB_1990.nc'];

varname = 'x';
x = ncread(source,varname)/1000; % easting-OSGB36 Grid reference
varname = 'y';
y = ncread(source,varname)/1000; % northing-OSGB36 Grid reference

dx = 1;
E = min(x):dx:max(x);
N = min(y):dx:max(y);
[X_coor,Y_coor] = meshgrid(E,N);
RAIN = [];

for year = YEARVEC
    for mon = 1:12

        source = ['K:\GEAR-month\CEH_GEAR_monthly_GB_',num2str(year),'.nc'];
        varname = 'rainfall_amount';
        RAIN_Aux=ncread(source,varname,...
            [1,1,mon],...
            [Inf,Inf,1]);
        RAIN_Aux = RAIN_Aux(:,end:-1:1,:);
        varname = 'min_dist';
        DIST=ncread(source,varname,...
            [1,1,mon],...
            [Inf,Inf,1]);
        
        
        RAIN_Aux = permute(RAIN_Aux,[2,1]);%[N,E]
        DIST = permute(DIST,[2,1]);
        RAIN(year-YEARVEC(1)+1,mon,:,:) = RAIN_Aux;

        
    end
end

if strcmpi(scale,'year')
    GEAR = squeeze(nansum(RAIN,2));
elseif strcmpi(scale,'month')
    GEAR = squeeze(nanmean(RAIN,1));
else
end

end