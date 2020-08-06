
region = getfield(REGIONS_info(),'WAL');
[HadUKMonth,E,N,HadUKAnnual] = getHadUK('tas');

[E1,N1] = getEN(region);
tas_JJA = squeeze(sum(HadUKMonth(:,:,6:8).*reshape([30,31,31],[1,1,3]),3))./92;
ii = E<=max(E1(:)) & E>=min(E1(:)) & N<=max(N1(:)) & N>=min(N1(:));
nanmean(tas_JJA(ii))

