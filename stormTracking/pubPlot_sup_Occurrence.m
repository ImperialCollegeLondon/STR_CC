% 
% get: Table of rainfall occurrence for each sub-region
% (data: is statistics computed for land areas)
% 
% Yuting Chen
% Imperial College London
% update: 2020/08/06

range = [0.1,10,30,50,Inf];
ENSEMBLENO = getEnsNos();

% occurrence of rainfall hourly events
range = [0.1,10,30,50,70,Inf];

tab = [];
rediff = [];
[tab(1:2,:),rediff(1,:)] = getTablePmax('CPM_NW',range);
[tab(3:4,:),rediff(2,:)] = getTablePmax('CPM_NE',range);
[tab(5:6,:),rediff(3,:)] = getTablePmax('CPM_S',range);

figure;hold on
ax = gca;
plot(rediff(1,:),'.-')
plot(rediff(2,:),'.-')
plot(rediff(3,:),'.-')
plot([ax.XLim],[0,0],'k--','linewidth',0.5);
setFigureProperty('Single')

ax.XTick = 1:length(range)-1;
ax.XTickLabels = {'0.1-10','10-30','30-50','50-80','>80'};
xtickangle(90)
ax.YLim = [-0.5,4];
ax.YTick = [-0.5,0,1,2,3,4];
ax.YTickLabels = {'-50%','0','100%','200%','300%','400%'};
ylabel('Relative Change of PDF')
xlabel('P_{max}(mm/h)')
legend('CPM-NW','CPM-NE','CPM-S')
set(gca,'linewidth',2)
% box on


%% AUXILLARY

function [tab,rediff] = getTablePmax(regionName,range)
[occ_past,occ_future] = deal([]);
ENSEMBLENO = getEnsNos();
for ensNo = ENSEMBLENO
    load(['D:\UKCP18\ConvectiveStorms\CS_',regionName,'_1980-2000_STATS_06-08_',ensNo{1},'.mat'])
    occ_past(end+1,:) = getPmax(STATS.rpmax,range);
    load(['D:\UKCP18\ConvectiveStorms\CS_',regionName,'_2060-2080_STATS_06-08_',ensNo{1},'.mat'])
    occ_future(end+1,:) = getPmax(STATS.rpmax,range);
end
tab = [nanmedian(occ_past,1);nanmedian(occ_future,1)];
tab(:,end+1) = nansum(tab(:,1:end),2);
rediff = nanmedian((occ_future-occ_past)./occ_past,1);
tab = tab/20;
end


function occ = getPmax(Pmax,range)

X = Pmax;
for ri = 1:length(range)-1
    occ(ri) = numel(X(X<range(ri+1) & X>=range(ri)));
end
end
