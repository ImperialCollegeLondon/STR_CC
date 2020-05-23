function [GEAR,X_coor,Y_coor,DIST] = getRegionalGEAR(region,scale,YEARVEC);
% GETREGIONALGEAR(region,scale) gived an average matrix, for a given scale,
%    given yearvec;
%
% Output format: GEAR 3d # [mon,northing,easting]
%
% Example: 
%         [GEAR,X_coor,Y_coor,DIST] = getRegionalGear(region,'month');
%
% @ yt


arguments
    region (1,1) struct
    scale (1,:) char = 'month'
    YEARVEC (1,:) double = 1980:1999;
end

dx = 1;
E = region.minE:dx:region.minE+(region.dimE-1)*region.dx;
N = region.minN:dx:region.minN+(region.dimN-1)*region.dx;

[X_coor,Y_coor] = meshgrid(E,N);

RAIN = NaN(0,0,size(X_coor,1),size(X_coor,2));
for year = YEARVEC
    for mon = 1:12
        [XX,YY,RAIN(year-YEARVEC(1)+1,mon,:,:),DIST] = import_GEAR_MONTHLY(X_coor*1000,Y_coor*1000,year,mon,false);
    end
end

GEAR = squeeze(nanmean(RAIN,1));

end