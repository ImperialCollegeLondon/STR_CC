
function [CPM,E,N] = getMonthCPMEns(region,mon)
% GETMONTHCPMENS(region,mon) give 12 ensembles mean pr for month <mon>
% Output format: CPM: {ensNo}([locE,locN,yrNo])
% <region> is not specified yet.
%


CPM_saved = load('D:\UKCP18\UK\MonMean_CPM.mat','RainEnsembles','E','N',...
    'EnsNo','IntFac','region');
OneMonth = CPM_saved.RainEnsembles{mon};
CPM = transpose(cellfun(@(OneEns)squeeze(nanmean(OneEns,3)),...
    OneMonth,'UniformOutput', false));% 12 ensembles

% CPM = cellfun(@(OneMon)transpose(cellfun(@(OneEns)...
%     squeeze(nanmean(OneEns,3)),...
%     OneMon,'UniformOutput', false)), ...
%     CPM_saved.RainEnsembles, 'UniformOutput', false);% get all month, all
%     ensembles
% CPM = cat(1,CPM{:});% to good cell format: [rowNo: mon; colNo: ensNo]

E = CPM_saved.E;
N = CPM_saved.N;

end
