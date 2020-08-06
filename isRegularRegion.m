function ret = isRegularRegion(region)
if isfield(region,'minE')
    ret = 1;
else
    ret = 0;
end
end