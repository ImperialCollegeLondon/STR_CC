
function [RAD,E,N] = getMonthRadar(region)
%
% output format: 
%               RAD: [mon, E, N];
%               E:
%               N:
% @yuting


dx = 2.2;
E = region.minE:dx:region.minE+(region.dimE-1)*region.dx;
N = region.minN:dx:region.minN+(region.dimN-1)*region.dx;
[N,E] = meshgrid(N,E);% not same as X-coor,Y-coor

SPath = 'D:/UKCP18/';
UK = load(sprintf('%sNIMRODRainPattern.mat',SPath),'MERain','E','N');%[Mon,Year,2D Matrix]

RAD = [];
in = getTrimTag('unit','km','product','radar2.2');

%[mon,enNo,[2d matrix]]
for mon = 1:12
    RAD(mon,:,:) = squeeze(nanmean(UK.MERain(mon,:,:,:),2));
    % mmrain{mon}(~in) = NaN;
end

indE = findClose(region.minE,UK.E(:,1));
indN = findClose(region.minN,UK.N(:,1));

RAD = RAD(:,...
    indE:indE+region.dimE-1,...
    indN:indN+region.dimN-1);

    function indX = findClose(x,Xvec)
        dis = abs(Xvec - x);
        indX = find(dis == min(dis(:)),1);
    end

end
