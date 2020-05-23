UKMap = getUKMap();
i0 = 1;
i1 = 1;
figure;hold on;
for i = 1:length(UKMap.borderE)
% plot([UKMap.borderE]/1000,[UKMap.borderN]/1000,'k');
if isnan(UKMap.borderE(i))
    i1 = i;
    fill([UKMap.borderE([i0:i1-1,i0])]/1000,[UKMap.borderN([i0:i1-1,i0])]/1000,ones(1,3)*0.8,'EdgeColor','none');
   
    i0 = i+1;
end

end
axis off
axis equal


% get several regions of interest
REGIONS = REGIONS_info();

plotOneRe = @(region)1;
% rectangle('position',[region.minE,region.minN,...
%     region.dimE*region.dx,...
%     region.dimN*region.dx],...
%     'linewidth',1,'Linestyle','-','EdgeColor',ones(1,3)*0.7);
tagOneRe = @(region)text(region.minE+region.dimE*region.dx/2,...
    region.minN+region.dimN*region.dx/2,region.Name,...
    'fontsize',8,'horizontalalignment','center','fontweight','bold',...
    'background',[0.95 0.95 0.95 0.1],'Color',ones(1,3)*0.5);

plotOneRe(REGIONS.SCO);
plotOneRe(REGIONS.WAL);
plotOneRe(REGIONS.EUK);

tagOneRe(REGIONS.SCO);
tagOneRe(REGIONS.WAL);
tagOneRe(REGIONS.EUK);

hold off
axis equal
xlim([0,800]);