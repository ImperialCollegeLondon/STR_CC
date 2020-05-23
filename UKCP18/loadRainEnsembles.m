function [E,N,RainEnsembles,scaleF,region] = loadRainEnsembles(region,MONS,Period)
% [E,N,RainEnsembles,scaleF,region] = loadRainEnsembles(region,MONS)
% is to : Load Rain Ensembles from UKCP18 
% (12 ensembles for several given month)
%
% Input: Region: <struct> info of region
%        MONS: <vector> scalar or vector within a range of [1,12];
%        Period: <char> '1980-2000';'2020-2040';'2060-2080'
% Output: ...
%
% @ Yuting Chen
% yuting.chen17@imperial.ac.uk

[~,~,RainEns_aux,~,~] = loadRegionFile(region,MONS(1),Period);
[RainEnsembles] = deal(cell(size(RainEns_aux,1),1));
RainEnsembles = cellfun(@(X)zeros(region.dimE,region.dimN,0),RainEnsembles,...
    'UniformOutput', false);
for mon = MONS
    fprintf('mon %2d\n',mon)
    [E,N,RainEns_aux,scaleF,region] = loadRegionFile(region,mon,Period);
    RainEnsembles = cellfun(@(X,Y)cat(3,X,Y),RainEnsembles,RainEns_aux,...
        'UniformOutput', false);
end

    function [E,N,RainEnsembles,scaleF,region] = loadRegionFile(region,mon,Period) %#ok<STOUT>
        [E,N] = getEN(region);
        % Compute the GridAreaMatrix
        region.gridArea = diff([E(1,1)-1.1;E(:,1)+1.1])...
            *diff([N(1,:)-1.1,N(1,end)+1.1]);
        
        if strcmp(Period,'1980-2000')
            filetil = 'D:/UKCP18/';
        elseif strcmp(Period,'2020-2040')
            filetil = 'D:/UKCP18_Future/';
        elseif strcmp(Period,'2060-2080')
            filetil = 'D:/UKCP18_Future_2060_2080/';
        end
        fileDir = [ filetil,region.Name,'/Ensems',...
            '_mon',sprintf('%02d',mon)];
        load([fileDir,'.mat'],'RainEnsembles');
        scaleF = 32;
    end
end