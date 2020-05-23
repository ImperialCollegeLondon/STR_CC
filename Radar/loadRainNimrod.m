%
% read Extracted Region Data
%

function [E,N,RainNimrods,scaleF,region] = loadRainNimrod(region,MONS)
% LOADRAINNIMROD is to load the radar data for a certain region, given
% months.
% NIMROD region file is required (which can be computed from 
% 'Radar_extractRegion.m'
%
% Input: region <struct> which is a field of the output of REGION_info();
%        MONS <vector> double
% Output:E <matrix> [meter]
%        N <matrix> [meter]
%        RainNimrods <3d matrix> in a direction of [E,N,Time]
%        scaleF <double> [32 as default]
%        region <struct>
%
% Example:
% REGIONS = REGIONS_info();
% region = REGIONS.London;
% MONS = [5,7];
% % save data for each month seperately;
% [E,N,RainNimrod,scaleF,region] = loadRainNimrod(region,MONS);
%
% @ Yuting Chen
% Imperial College London
% Update: 2019.12.06


RainNimrods = int16(NaN(region.dimE,region.dimN,0));

for mon = MONS
    fprintf('mon %2d\n',mon)
    [E,N,Rain_aux,scaleF,region] = loadRegionFile(region,mon);
    RainNimrods = cat(3,RainNimrods,Rain_aux);
end


% Compute the GridAreaMatrix
region.gridArea = 2.2.^2*ones(size(RainNimrods,1),size(RainNimrods,2));


    function [E,N,Rain,scaleF,region] = loadRegionFile(region,mon) %#ok<STOUT>
        
        [E,N] = getEN(region);
        
        fileDir = ['D:/UKCP18/',region.Name,'/Radar',...
            '_mon',sprintf('%02d',mon)];
        load([fileDir,'.mat'],'Rain');
        scaleF = 32;
    end

end



% load('Radar_2007_London.mat')
% imagesc(squeeze(nanmean(PRS2_2,3))/scaleF*24*365);
% EE_del = imresize(RADAR_info.RegionEE,[50,50]);
% NN_del = imresize(RADAR_info.RegionNN,[50,50]);
% for i = 1:8760-8757
%     contourf(EE_del,NN_del,squeeze(PRS2_2(:,:,i))/scaleF,[0:01:20]);
%     colormap(pink)
%     axis off
%     box on
%     pause(0.2);
% end


