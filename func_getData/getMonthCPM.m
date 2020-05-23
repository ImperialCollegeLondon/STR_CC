
function [CPM,E,N] = getMonthCPM(region)
% GETMONTHCPM(region) gives an average mean monthly pr (average over 12 ens)
% Input:
% 
% Output:
% 
% Output format: {monNo,ensNo}[E,N]
%
% Example:
%        [CPM,E,N] = getMonthCPM(region);
% @ Yuting

arguments
    
    region (1,1) struct
    
end

if strcmpi(region.Name,'UK')
    
    CPM_saved = load('D:\UKCP18\UK\MonMean_CPM.mat','RainEnsembles','E','N',...
        'EnsNo','IntFac','region');
    CPM = cellfun(@(OneMon)transpose(cellfun(@(OneEns)...
        squeeze(nanmean(OneEns,3)),...
        OneMon,'UniformOutput', false)), ...
        CPM_saved.RainEnsembles, 'UniformOutput', false);
    CPM = cat(1,CPM{:});
    E = CPM_saved.E;
    N = CPM_saved.N;
    
else %extract local mat for <region>
    
    CPM_saved = load('D:\UKCP18\UK\MonMean_CPM.mat','RainEnsembles','E','N',...
        'EnsNo','IntFac','region');
    
    [region.i,region.j] = getRegionIJ(CPM_saved.E,CPM_saved.N,region.minE,region.minN);
    
    % E1 = E(region.i:region.i+region.dimE-1,region.j);
    % N1 = N(region.i,region.j:region.j+region.dimN-1);
    
    
    CPM = cellfun(@(OneMon)transpose(cellfun(@(OneEns)...
        squeeze(nanmean(...
        squeeze(OneEns(region.i:region.i+region.dimE-1,...
        region.j:region.j+region.dimN-1,:)),...
        3)),...
        OneMon,'UniformOutput', false)), ...
        CPM_saved.RainEnsembles, 'UniformOutput', false);
    CPM = cat(1,CPM{:});
    E = CPM_saved.E(region.i:region.i+region.dimE-1,region.j);
    N = CPM_saved.N(region.i,region.j:region.j+region.dimN-1);
    
    [N,E] = meshgrid(N,E);
    
end

end



% - Below: Save in a cell format -- Time CONSUMING
% REGIONS = REGIONS_info();
% region = REGIONS.UK;
% 
% EnsNo=getEnsNos();
% % EnsNo = EnsNo(1:2);
% IntFac = 1;
% RainEnsembles = [];
% 
% for MON = 1:12
%     tic
%     [RainEnsembles{MON},IntFac,E,N] = readCPM_pr_average(region,EnsNo,MON,IntFac);
%     toc
%     fprintf('Mon%02d done',MON);
% end
% 
% save('D:\UKCP18\UK\MonMean_CPM.mat','RainEnsembles','E','N',...
%     'EnsNo','IntFac','region','-v7.3'); % Space needed: 5.97GB
%
%
% -- Below: Save in a vector format -- Time CONSUMING
%
% CPM = cellfun(@(OneMon)cell2mat(cellfun(@(OneEns)...
%     reshape(OneEns,[1,size(OneEns)]),...
%     OneMon,'UniformOutput', false)), ...
%     CPM_saved.RainEnsembles, 'UniformOutput', false);% each cell represent one month
% 
% CPM = cell2mat(cellfun(@(OneMon)reshape(OneMon,[1,size(OneMon)]), ...
%     transpose(CPM), 'UniformOutput', false));
% 
% save('D:\UKCP18\UK\Mon_Ens_E_M_Yr_CPM.mat','CPM','E','N',...
%     'EnsNo','IntFac','region','-v7.3')
%
