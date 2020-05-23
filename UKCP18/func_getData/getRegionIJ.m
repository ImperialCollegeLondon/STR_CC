function [i,j] = getRegionIJ(E,N,e0,n0)
% for a format: [E,N]
%
if isvector(E)
    [N,E]=meshgrid(N,E);
end

disKm = (E-e0).^2+(N-n0).^2;
[i,j] = find(disKm == min(disKm(:)),1);

if ~(isscalar(i) & isscalar(j))
    error('Check i,j');
end

end