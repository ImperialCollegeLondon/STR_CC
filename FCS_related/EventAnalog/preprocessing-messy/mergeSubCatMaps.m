% ----------------------------------------------------------------------- %
% Get Flood Maps for each subwater course WC$**$
%
% Squeeze image to vector
% Yuting Chen
% Update: 2020.02.01
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
%%
filename = [sPath,'Birmingham+RPS_WC'];
savePlot(filename,'XYWH',[150,0,400,400]);%,'needreply','N');


%% Get and Merge Info of FloodMap Files for this EventNo

D = dir([FloodPath,'Batch_*\*_WC*_PEAK.asc']);

EventNames = unique(cellstr(arrayfun(@(x)x.name(1:3), ...
    D, 'UniformOutput', false)));

WCNames = unique(cellstr(arrayfun(@(x)regexp(string(x.name(7:8)),'\d*','match'), ...
    D, 'UniformOutput', false)));

for WCNo = WCNames'
    
    WCNo = WCNo{1};
    FloodFNs = getWCFiles(FloodPath,WCNo);
    
    % Merge RPS Files
    
    [FVec,FMap,FMapFName,FMapInfo,E_rps,N_rps] = SqueezeeRPSFloodMap(FloodFNs,WCNo);
    
    unit = 'meter';
    
end

%% add those event without WC info as 0;

for WCNo = WCNames'
    WCNo = WCNo{1};
    for EventNo = EventNames'
        
        EventNo = EventNo{1};
        filename = sprintf('G:/BIGDATA/TOPIC 2/ProcessedFiles/WCMaps/WCMaps_Event%s_WCNo%02d.mat',...
            EventNo,str2num(WCNo));
        %'FloodVec','CC','EVec','NVec','unit');
        fid = fopen(filename);
        if fid == -1
            FNs = dir(sprintf('G:/BIGDATA/TOPIC 2/ProcessedFiles/WCMaps/WCMaps_*_WCNo%02d.mat',...
                str2num(WCNo)));
            filename = [FNs(1).folder,filesep,FNs(1).name];
            A = load(filename);
            FloodVec = A.FloodVec*0;
            CC = A.CC;
            EVec = A.EVec;
            NVec = A.NVec;
            unit = A.unit;
            save(filename,'FloodVec','CC','EVec','NVec','unit')
        else
            fclose(fid);
        end
        
    end
end




% [row_nort,col_east] = find(FMap);
% subplot 121
% plot(col_east,row_nort,'.');
% subplot 122
% plot(sparse(FMap));


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


function FileNames = getWCFiles(Path,TagNo)


D = dir([Path,'Batch_*\*_WC',TagNo,'_*PEAK.asc']);
FileNames = arrayfun(@(x)strcat(string(x.folder),'\',string(x.name)), ...
    D, 'UniformOutput', false);

end


function plotWCs(Path)

TagNo = '*';
D = dir([Path,'\Batch*\',TagNo,'_WC*_PEAK.asc']);
FileNames = arrayfun(@(x)strcat(string(x.folder),'\',string(x.name)), ...
    D, 'UniformOutput', false);

for fn = 1:numel(FileNames)
    thisRec = getRec(FileNames{fn});
    rectangle('Position',[thisRec.x(1) thisRec.y(1) thisRec.dx thisRec.dy],...
        'LineWidth',2,'EdgeColor','k')
    % text(thisRec.x(1)+500,thisRec.y(1)+1000,regexp(FileNames{fn},'WC\d*','match'));
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



function [FVec,FMap,FTag,FMapInfo,EH,NH] = SqueezeeRPSFloodMap(FloodFNs,WCNo)

% this will convert 

try
    
    Birm = getfield(REGIONS_info(),'Birm');
    % [E,N] = meshgrid(round(1000*(Birm.minE:0.01:Birm.minE+Birm.dx*(Birm.dimE-1))),...
    %     round(1000*flip(Birm.minN:0.01:Birm.minN+Birm.dx*(Birm.dimN-1))));
    % FMap = NaN(size(E)); % dim:[N]*[E]
    % [FTag,FMapInfo] = deal([]);
    [FMap,FTag,FMapInfo] = deal([]);
    
    for Num = 1:numel(FloodFNs)
        
        info = struct;
        
        [OUT, info.ncols, info.nrows, info.xllcorner, info.yllcorner, info.cellsize, info.nodata] = ...
            ascii_reader(FloodFNs{Num});
        
        [EH,NH] = meshgrid(info.xllcorner:info.cellsize:info.xllcorner+(info.ncols-1)*info.cellsize,...
            flip(info.yllcorner:info.cellsize:info.yllcorner+(info.nrows-1)*info.cellsize));

        FTag{Num} = FloodFNs{Num};
        FMapInfo{Num} = info;
        
        if size(FMap)~=size(OUT)
            FMap = zeros(size(OUT));
        else
            FMap = FMap+OUT;
        end
        
        fprintf('%s - analyse - done.\n',FloodFNs{Num})
    end
    
catch me
    1;
end


bw = FMap > 0;

CC = bwconncomp(bw);
Ind = cell2mat(CC.PixelIdxList');
EVec = EH(Ind);
NVec = NH(Ind);
FVec = [];

% start to squeeze

for Num = 1:numel(FloodFNs)
    
    info = struct;
    
    [OUT, info.ncols, info.nrows, info.xllcorner, info.yllcorner, info.cellsize, info.nodata] = ...
        ascii_reader(FloodFNs{Num});
    
    FVec(Num,:) = reshape(OUT(Ind),1,[]);
    
    fprintf('%s - squeeze - done.\n',FloodFNs{Num})
    
    EventNo = regexp(regexp(FloodFNs{Num, 1},'\d*_WC','match'),'\d*','match');
    FloodVec = FVec(Num,:)';
    unit = 'meter';
    
    save(sprintf('G:/BIGDATA/TOPIC 2/ProcessedFiles/WCMaps/WCMaps_Event%s_WCNo%02d.mat',...
        EventNo,str2num(WCNo)),'FloodVec','CC','EVec','NVec','unit');
    
end

end









