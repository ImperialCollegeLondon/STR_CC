function in = getTrimTag(NameValueArgs);
% GETRIMTAG will return the corresponding logical tag to trim data.
% Example: 
%        in = getTrimTag('unit','km','product','radar2.2');
%        
%        in = getTrimTag('unit','km','product','cpm2.2');
% @ Yuting
% Imperial College London

arguments
    NameValueArgs.unit
    NameValueArgs.product
end

% inpolygon(E,N,UKMap.borderE/1000,UKMap.borderN/1000);

if strcmpi(NameValueArgs.unit,'KM') && strcmpi(NameValueArgs.product,'RADAR2.2')
    
    load('H:\CODE_MATLAB\SpatialTemporalDATA\shapeFileFolder\radar_trimTag.mat','in');
    
elseif strcmpi(NameValueArgs.unit,'KM') && strcmpi(NameValueArgs.product,'cpm2.2')

    load('H:\CODE_MATLAB\SpatialTemporalDATA\shapeFileFolder\cpm_trimTag.mat','in');
    
else
    
    error('havent presave trimMat for other product.')
end


end