function [E1,N1] = getEN(region)

filePath = ['K:/UkCp18/01/pr_rcp85_land-cpm_uk_2.2km_'];
ncFileName = [filePath,'01','_1hr_19910801-19910830.nc'];
LAT=ncread(ncFileName,'latitude');
LON=ncread(ncFileName,'longitude');

[E, N] = ll2os(LAT, LON);
E = E/1000;
N = N/1000;
if isRegularRegion(region)
    [region.i,region.j] = getRegionIJ(E,N,region.minE,region.minN);
    E1 = E(region.i:region.i+region.dimE-1,region.j);
    N1 = N(region.i,region.j:region.j+region.dimN-1);
    [E1,N1] = meshgrid(E1,N1);
    E1 = E1';
    N1 = N1';
else
    angle = atan((region.E(2)-region.E(1))/(region.N(2)-region.N(1)))*180/pi;
    E = imrotate(E,angle,'nearest','loose');
    N = imrotate(N,angle,'nearest','loose');
    [region.i,region.j] = arrayfun(@(e0,n0)getRegionIJ(E,N,e0,n0),region.E,region.N);
    i_this = min(region.i):max(region.i);
    j_this = min(region.j):max(region.j);
    E1 = E(i_this,j_this);
    N1 = N(i_this,j_this);
end
end
