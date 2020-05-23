% ----------------------------------------------------------------------- %
% Merge Flood maps for each rainfall event
%
% Convert to a better storing format
% .nc or .mat or single/double/int<*10?>
%
% Yuting Chen
% Update: 2020.01.31
% ----------------------------------------------------------------------- %


%% Configuration

FloodPath = 'K:\DATA_FCS\RPS_FloodMaps\';
RainfallPath = 'G:\BIGDATA\TOPIC 2\FileFromLIPEN\Rainfall_in_157Floods\KED_ASCII\';


EventNo = '037';
EventTimes = '201408011300-201408021535';



%% Get Info for Birmingham District

fileName = 'K:\UK_shape\bdline_essh_gb\Data\GB\district_borough_unitary_region.shp';
area = shaperead(fileName);%,'Attributes',{'BoundingBox','X','Y','Name'});
area = area(109); % name: {'Birmingham District (B)'}
REGIONS = REGIONS_info();
region = REGIONS.Birm;
[Topo,E,N] = getTopo(region);


%% plot Info

pcolor(E,N,Topo);shading flat
cptcmap('DEM_print', 'mapping','direct');
caxis([0,500])

hold on;
plotWCs(FloodPath)
c = area.Y;
plot(area.X,area.Y,'r-','linewidth',4);
colorbar

filename = [sPath,'Birmingham+RPS_WC'];
savePlot(filename,'XYWH',[150,0,400,400]);%,'needreply','N');


%% Get and Merge Info of FloodMap Files for this EventNo

D = dir([FloodPath,'\Batch*\*_WC*_PEAK.asc']);

EventNames = unique(cellstr(arrayfun(@(x)x.name(1:3), ...
    D, 'UniformOutput', false)));

for EventNo = EventNames'
    
    EventNo = EventNo{1};
    FloodFNs = getFloodFiles(FloodPath,EventNo);
    
    % Merge RPS Files
    
    [FMap,FMapFName,FMapInfo,E,N] = ImportRPSFloodMap(FloodFNs);
    
    unit = 'cm';
    
    save(sprintf('G:/BIGDATA/TOPIC 2/ProcessedFiles/FloodMaps_40_30/FloodMaps_Merged_EventNo%s.mat',...
        EventNo),'FMap','FMapFName','FMapInfo','E','N','unit');
    
end



% [row_nort,col_east] = find(FMap);
% subplot 121
% plot(col_east,row_nort,'.');
% subplot 122
% plot(sparse(FMap));

% add those event without flood
for EventNo = 1:157
    
    try
        load(sprintf('G:/BIGDATA/TOPIC 2/ProcessedFiles/FloodMaps_40_30/FloodMaps_Merged_EventNo%03d.mat',...
            EventNo),'FMap','FMapFName','FMapInfo','E','N','unit');
    catch me 
        FMap = 0*E;
        FMapFName = NaN;
        FMapInfo = NaN;
        save(sprintf('G:/BIGDATA/TOPIC 2/ProcessedFiles/FloodMaps_40_30/FloodMaps_Merged_EventNo%03d.mat',...
            EventNo),'FMap','FMapFName','FMapInfo','E','N','unit');
    end
    
end

 
%% extract the flood maps having same location+same size of rain maps
for EventNo = 1:157
    try

        load(sprintf('G:/BIGDATA/TOPIC 2/ProcessedFiles/FloodMaps_40_30/FloodMaps_Merged_EventNo%03d.mat',...
            EventNo),'FMap','FMapFName','FMapInfo','E','N','unit');
        fprintf('%03d: Done \n',EventNo);
        
        FMap = FMap(end-33:end,:);
        E = E(end-33:end,:);
        N = N(end-33:end,:);
        N = flip(N,1);
        FMap = flip(FMap,1);
        save(sprintf('G:/BIGDATA/TOPIC 2/ProcessedFiles/FloodMaps_34_30/FloodMaps_Merged_EventNo%03d.mat',...
            EventNo),'FMap','FMapFName','FMapInfo','E','N','unit');
        
        
    catch me
        fprintf('No Event Info this time \n');
    end
    
end



%% Show FloodMap

INDALL = [];

for eventNo = 1:157
    
    EVENT = sprintf('%03d',eventNo);
    
    try
    load(sprintf('G:/BIGDATA/TOPIC 2/ProcessedFiles/FloodMaps/FloodMaps_Merged_EventNo%s.mat',...
        EVENT),'FMap','FMapFName','FMapInfo','E','N','unit');
    
    
    hazardR = FMap/100*(0.25+0.5)+...
        (FMap/100<0.25)*0.5+(FMap/100>=0.25)*1;

    
    [i,j] = find(FMap>10);
    IND = sub2ind(size(E),i,j);
    scatter(E(IND),N(IND),20,'r','fill');
    alpha(0.2)
    drawnow
    hold on;
    
    INDALL = [INDALL;IND];
    
    catch me
        fprintf('No event %03d\n',eventNo);
    end
    
end
%%
[Nind,Eind] = ind2sub(size(E),INDALL);
[Block,XEDGES,YEDGES] = histcounts2(Nind,Eind,0.5:size(N,1)+0.5,0.5:size(E,2)+0.5);
% pcolor(exp(FMap)-1);
% shading flat;
% cptcmap('GMT_drywet', 'mapping', 'scaled');
% colorbar;
% pause(0.5);%

f = @(x)imresize(x,0.02,'box')*50;
pcolor(f(E),f(N),f(Block));
shading flat
cptcmap('GMT_haxby', 'mapping','scaled');

%% AUXILLARY FUNCTION


function [Topo,E,N] = getTopo(region)

TERRAIN = load(['K:\UK_shape\DTM50.mat'],'DTM50','Eno','Nno');

Nind_max = findClose(region.minN*1000,TERRAIN.Nno);
Nind_min = findClose(1000*(region.minN+region.dx*(region.dimN-1)),TERRAIN.Nno);

Eind_min = findClose(region.minE*1000,TERRAIN.Eno);
Eind_max = findClose(1000*(region.minE+region.dx*(region.dimE-1)),TERRAIN.Eno);

Topo = TERRAIN.DTM50(Nind_max:-1:Nind_min,...
    Eind_min:Eind_max);

[NNno,EEno] = meshgrid(TERRAIN.Nno(Nind_max:-1:Nind_min),TERRAIN.Eno(Eind_min:Eind_max));
E = EEno;
N = NNno;
Topo = transpose(Topo);

    function indX = findClose(x,Xvec)
        dis = abs(Xvec - x);
        indX = find(dis == min(dis(:)),1);
    end

end


function FileNames = getFloodFiles(Path,TagNo)


D = dir([Path,'\Batch*\',TagNo,'_WC*_PEAK.asc']);
FileNames = arrayfun(@(x)strcat(string(x.folder),'\',string(x.name)), ...
    D, 'UniformOutput', false);

end


function plotWCs(Path)

TagNo = '*';
D = dir([Path,TagNo,'_WC*_PEAK.asc']);
FileNames = arrayfun(@(x)strcat(string(x.folder),'\',string(x.name)), ...
    D, 'UniformOutput', false);

for fn = 1:numel(FileNames)
    thisRec = getRec(FileNames{fn});
    rectangle('Position',[thisRec.x(1) thisRec.y(1) thisRec.dx thisRec.dy],...
        'LineWidth',2,'EdgeColor','k')
    hold on;
end

xlabel('Easting')
ylabel('Northing')

    function thisRec = getRec(filename)
        
        fin = fopen(filename,'r');
        
        A = fscanf(fin,'%s',1); ncols = fscanf(fin,'%f',1);           %#ok<NASGU>
        A = fscanf(fin,'%s',1); nrows = fscanf(fin,'%f',1);           %#ok<NASGU>
        A = fscanf(fin,'%s',1); xllcorner = fscanf(fin,'%f',1);       %#ok<NASGU>
        A = fscanf(fin,'%s',1); yllcorner = fscanf(fin,'%f',1);       %#ok<NASGU>
        A = fscanf(fin,'%s',1); cellsize = fscanf(fin,'%f',1);        %#ok<NASGU>
        
        thisRec = struct;
        
        thisRec.x = [xllcorner,xllcorner+(nrows-1)*cellsize];
        thisRec.y = [yllcorner,yllcorner+(ncols-1)*cellsize];
        thisRec.d0 = cellsize;
        thisRec.dx = (nrows-1)*cellsize;
        thisRec.dy = (ncols-1)*cellsize;
    end

end



function [FMap,FTag,FMapInfo,E,N] = ImportRPSFloodMap(FloodFNs)

try
    
    Birm = getfield(REGIONS_info(),'Birm');
    reso = 1000;%1000 m
    [E,N] = meshgrid(round(1000*(Birm.minE:reso/1000:Birm.minE+Birm.dx*(Birm.dimE-1))),...
        round(1000*flip(Birm.minN:reso/1000:Birm.minN+Birm.dx*(Birm.dimN-1))));
    FMap = zeros(size(E)); % dim:[N]*[E]
    [FTag,FMapInfo] = deal([]);
    
    for Num = 1:numel(FloodFNs)
        
        info = struct;
        
        [OUT, info.ncols, info.nrows, info.xllcorner, info.yllcorner, info.cellsize, info.nodata] = ...
            ascii_reader(FloodFNs{Num});
        
        [OUT] = aggregateImage(OUT*100,round(1000/info.cellsize),'max');
        
        [info] = updateInfo(size(OUT),info,1000);
        
        % FMap{Num} = imresize(OUT,1/25);
        
        [r1,r2,c1,c2] = findLoc(info,E,N);
        
        FMap(r1:r2,c1:c2) = OUT;
        
        FTag{Num} = FloodFNs{Num};
        FMapInfo{Num} = info;
        
        fprintf('%s done.\n',FloodFNs{Num})
    end
    
catch me
    1;
end
    function [r1,r2,c1,c2] = findLoc(info,E,N)
        [r2,c1] = find(E == info.xllcorner & N == info.yllcorner);
        r1 = r2-(info.nrows-1);
        c2 = c1+(info.ncols-1);
    end

    function [info] = updateInfo(matSize,info,newCellSize)
        info.ncols = matSize(2);
        info.nrows = matSize(1);
        info.cellsize = newCellSize;
        info.xllcorner = round(info.xllcorner/newCellSize)*newCellSize;
        info.yllcorner = round(info.yllcorner/newCellSize)*newCellSize;
    end

end









