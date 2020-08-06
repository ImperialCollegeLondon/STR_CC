function Area = getSumArea(regionName,in)
% unit: [Km^2]
% suitable for UKCP CPM-12 data
arguments
regionName (1,:) char
in (:,:) = ''
end
region = getfield(REGIONS_info(),regionName);
[E1,N1] = getEN(region);
Area = numel(E1(:))*2.2*2.2;
if ~isempty(in)
Area = numel(E1(in))*2.2*2.2;
end
end