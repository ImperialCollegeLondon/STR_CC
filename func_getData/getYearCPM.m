
function [CPM_year,E,N] = getYearCPM(region)
% GETYEARCPM(region)
% order:
% CPM_year: [Ens,E,N,yr]
%

CPM_saved = load('D:\UKCP18\UK\Mon_Ens_E_M_Yr_CPM.mat','CPM','E','N',...
    'EnsNo','IntFac','region');

E = CPM_saved.E;
N = CPM_saved.N;


[yy,mm] = meshgrid(1980:1999,1:12);
days = reshape(eomday(yy,mm),[12,1,1,1,20]);

CPM_year = CPM_saved.CPM*24.*days;
clear CPM_saved


CPM_year = squeeze(nansum(CPM_year,1));



% % OneMonth = CPM_saved.RainEnsembles{mon};
% % CPM = transpose(cellfun(@(OneEns)squeeze(nanmean(OneEns,3)),...
% %     OneMonth,'UniformOutput', false));% 12 ensembles



end
