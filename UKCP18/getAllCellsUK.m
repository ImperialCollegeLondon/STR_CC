function [E,N,CELLALL,THRE] = getAllCellsUK(season,region,source,cellThre)
% GETALLCELLSUK is to load saved cell info.
% source: 'CPM' or 'RAD'
% example:
%       [E,N,CELLALL] = getAllCellsUK(season,region.Name,source);
% @ Yuting Chen
% Imperial College London


arguments
    
    season (1,1) double
    region (1,1) struct
    source (1,:) char {mustBeMember(source,{'CPM','RAD'})} = 'CPM'
    cellThre (1,1) double = 1
    
end

[E,N] = getEN(source,region);
switch(source)
    
    case 'CPM'
        dataSP = 'H:\DATA_CLIMATE\UKCP18\';
        
        
        CELLALL = table;
        ENSEMBLENO = getEnsNos();
        for enNo = 1:length(ENSEMBLENO)
            for mon = getMons(season)
                cellFileName = sprintf('%sCellProp_UKCP%s_Mon%01d_%s.mat',...
                    dataSP,ENSEMBLENO{enNo},mon,region.Name);
                load(cellFileName,'CELLA','THRE');
                tab = CELLA{cellThre};
                tab.ensNo = repmat(enNo,size(tab,1),1);
                CELLALL = [CELLALL;tab];
            end
        end
    case 'RAD'
        dataSP = 'H:\DATA_CLIMATE\UKCP18\';
        CELLALL = table;
        for mon = getMons(season)
            cellFileName = sprintf('%sCellProp_RAD_Mon%01d_%s.mat',...
                dataSP,mon,region.Name);
            load(cellFileName,'CELLA','THRE');
            CELLALL = [CELLALL;CELLA{cellThre}];
        end
end

THRE = THRE(cellThre);


end



function [E,N] = getEN(source,region)
switch(source)
    case 'CPM'
        fileName = sprintf('pr_rcp85_land-cpm_uk_2.2km_%s_1hr_%04d%02d01-%04d%02d30.nc',...
            '01',1981,1,1981,1);
        filePath = ['K:/UkCp18/','01','/'];
        listaRain = [filePath,fileName];
        A = ncinfo(listaRain);
        
        LAT = ncread(listaRain,'latitude');
        LON = ncread(listaRain,'longitude');
        [E, N] = ll2os(LAT, LON);
        E = E/1000;
        N = N/1000;
    case 'RAD'
        E = region.minE:2.2:region.minE+(region.dimE-1)*2.2;
        N = region.minN:2.2:region.minN+(region.dimN-1)*2.2;
        [N,E] = meshgrid(N,E);
end
end



