
function [sPath] = StartUp()

addpath(genpath('C:\Users\Yuting Chen\Dropbox (Personal)\MatLab_Func'));
addpath(genpath(cd));
addpath(genpath(('H:\CODE_MATLAB\SpatialTemporalDATA\shapeFileFolder')));

Birm = struct;

% IMPORT: info from https://gridreferencefinder.com/
Birm.Easting = 406689;
Birm.Northing = 286822;

save([cd,'\Birm.mat'],'Birm');

sPath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_FCS\';
han = setFigureProperty('Meeting');

end

