% -------------------------------------------------------------------------- %
% THIS FILE IS TO PLOT DEMO RESULT WITH BACKGROUND OBTAINED FROM OS LOCAL MAP
% FLOOD DEPTHS (PEAK) ARE USED TO REPRESENT FLOOD RISK.
% SHOWN AS SCATTER POINTS
% SIZE/COLORMAP: REPRESENT THE MAGNITUDE OF DEPTH
% @ YUTING CHEN
% IMPERIAL COLLEGE LONDON
% yuting.chen17@Imperial.ac.uk
% -------------------------------------------------------------------------- %


% LOAD DATA
expNo = 'C1P5In-Rad1';eventNo4val = 90;evtPlot = 70;
eviPlot = eventNo4val;
[version,testConfig,filefolder] = getExpNoInfo(expNo);
load([filefolder,sprintf('STATS_%s_ev%03d',version,eventNo4val)],'STATS','FOREOUTPUT');
floodmaps_test = squeeze(FOREOUTPUT.FLO_test.floodmaps(evtPlot,:,:));
floodmaps_pred = squeeze(FOREOUTPUT.FLO_pred.floodmaps(evtPlot,:,:,:));
E = FOREOUTPUT.E;
N = FOREOUTPUT.N;

% PLOT BACKGROUND
[ax1] = plotBackground();

% PLOT FLOODMAP
floodobs = floodmaps_test;
floodpreds = floodmaps_pred;
axislim = nanmax(floodobs(:));

floodval = prctile(floodpreds,90,3);% floodobs;%
plotFlodeoth(ax1,axislim,floodval,E,N);

savePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_FCS';
fileTag = sprintf('%s_%s_evi%03d_ti%03d',version,testConfig.imageMethod,eviPlot,evtPlot);
filename = [savePath,filesep,'Prediction_DemoRes_FloodMap_',fileTag,'_PredPrct90'];
savePlot(filename,'XYWH',[50,50,400,300],'needreply','Y');


%% AUXILLARY FUNCTION

function [] = plotFlodeoth(ax1,axislim,floodval,E,N)

view(2)
ax2 = axes;
aggScale = 3; mode = 'max';
maps = aggregateImage(floodval,aggScale,mode);
E0 = imresize(E,'scale',1/aggScale,'method','box');
N0 = imresize(N,'scale',1/aggScale,'method','box');
plot_DemoFloodMaps(maps,E0,N0,ax1,ax2,axislim)

    function plot_DemoFloodMaps(maps,E0,N0,ax1,ax2,axislim)
        maps(maps == 0) = NaN;
        % pcolor(ax2,E,N,maps);shading flat;
        
        [cmap, lims, ticks, bfncol, ctable] = cptcmap('flood_blue','mapping','scaled','ncol',15);
        colfmap = NaN(numel(maps(:)),3);
        colfmap(~isnan(maps(:)),:) = cmap(getLevel(maps(~isnan(maps)),linspace(0,axislim,15)),:);
        scatter(ax2,E0(:),N0(:),(maps(:))*4,colfmap,'fill');shading flat;
        alpha(0.6)
        %%Link them together
        linkaxes([ax1,ax2])
        %%Hide the top axes
        ax2.Visible = 'off';
        ax2.XTick = [];
        ax2.YTick = [];
        
        %%Give each one its own colormap
        cptcmap('GMT_gray',ax1, 'mapping', 'scaled','flip',false);
        cptcmap('flood_blue', ax2,'mapping','scaled','ncol',15);
        axis off
        box on
        axis equal
        caxis(ax2,[0,axislim]);
        caxis(ax1,[80,255]);
        c = colorbar('location','Manual', 'position', [0.25 0.15 0.02 0.75]);
        t0 = c.TickLabels;
        c.TickLabels = strcat(t0,'cm');
    end

end

function [ax1] = plotBackground()
tic
f = figure;
f.GraphicsSmoothing = 'off';
ax1 = axes;
load('G:\BIGDATA\TOPIC 2\BirmLocalMap_reso100m.mat','x','y','imageData');
Area = getBoundaryShp('Birm');

XLIM = [Inf,0];YLIM = [Inf,0];

for i = 1:length(x)
    for tileNum = 1:numel(x{i})
        pcolor(x{i}{tileNum},y{i}{tileNum},imageData{i}{tileNum}); shading flat; hold on;
        XLIM = [min(XLIM(1),min(x{i}{tileNum})),max(XLIM(2),max(x{i}{tileNum}))];
        YLIM = [min(YLIM(1),min(y{i}{tileNum})),max(YLIM(2),max(y{i}{tileNum}))];
        
    end
    cptcmap('GMT_gray','mapping', 'scaled','flip',false,'ncol',256); caxis([0,255])
    drawnow;
end
plot(Area.X,Area.Y,'k-','linewidth',2)
set(gca,'YDir','normal')
axis equal
axis off
xlim(XLIM);ylim(YLIM)
toc
end

function mapLevel = getLevel(maps,lowThre)
mapLevel = NaN(size(maps));
lowThre = [lowThre,inf];
for li = 1:numel(lowThre)-1
    mapLevel(maps>=lowThre(li) & maps<lowThre(li+1)) = li;
end
end

function [Eran,Nran,Area] = loadBirmENRange()
Area = getBoundaryShp('Birm');
Eran = [min(Area.X),max(Area.X)];
Nran = [min(Area.Y),max(Area.Y)];
end

function [files] = getAllFiles(fn)
files = dir([fn,'*\data\*.tif']);
end

function [inTag] = isInBirm(Eran,Nran,FILES)

inTag = [];
for fi = 1:length(FILES)
    t = Tiff([FILES(fi).folder,'\',FILES(fi).name],'r');
    strc = imfinfo([FILES(fi).folder,'\',FILES(fi).name]);
    x = strc.ModelTiepointTag(4);
    y = strc.ModelTiepointTag(5);
    inTag(fi) = x<Eran(2)+5000&x>Eran(1)-5000&y<Nran(2)+5000&y>Nran(1)-5000;
end
inTag = logical(inTag);

end

function [reg,num] = extractFile(file)

fn = file.name;
reg = fn(1:2);
num = fn(3:4);

end

function [X,Y,imageDataS] = readOneTile(area,num,pl,reso)

filePath = 'K:\UK_shape\OC_UKOpenMap - Local\';
FILES = dir([filePath,'omlras_gtfc_',area,'\Data\',upper(area),num,'*.tif']);

[X,Y,imageDataS] = deal([]);
for fi = 1:numel(FILES)
    
    %%%% need to find out loaction of the tile.
    
    t = Tiff([FILES(fi).folder,'\',FILES(fi).name],'r');
    strc = imfinfo([FILES(fi).folder,'\',FILES(fi).name]);
    x = strc.ModelTiepointTag(4);
    y = strc.ModelTiepointTag(5);
    imageData = read(t);
    unit = reso;
    imageData = imresize(imageData,'scale',1/(unit-0.01),'method','box');
    
    x = (x+unit/2):unit:(x+unit*(size(imageData,1)-1)+unit/2);
    y = (y-unit/2):-unit:(y-unit*(size(imageData,1)-1)-unit/2);
    
    if pl
        % pcolor(x,y,imageData); shading flat;
        % cptcmap('GMT_gray','mapping', 'scaled','flip',false,'ncol',256); caxis([0,255])
        imshow(imageData,'XData',x,'YData',y) %'Colormap', summer(round((rand*200))));
        hold on;
        colormap(gca, gray(256))
        drawnow
    else
    end
    
    X{fi} = x;
    Y{fi} = y;
    imageDataS{fi} = imageData;
    
end
% axis normal

end




% ---------------------------------------------------------------------- %
% % Find out Birm map
% % WHICH HAS BEEN SAVED
%
% [Eran,Nran] = loadBirmENRange();
% [files] = getAllFiles('K:\UK_shape\OC_UKOpenMap - Local\omlras_gtfc_');
% [inTag] = isInBirm(Eran,Nran,files);
% files = files(logical(inTag));
% save('BirmLocalMapFiles.mat','files');

% % USING SAVED MAP TO PLOT BACKGROUND
% % WHICH HAS BEEN SAVED
%
% load('BirmLocalMapFiles.mat','files');
% reso = 100;
% pl = false;
% [x,y,imageData] = deal([]);
% for i = 1:length(files)
%     [reg,num] = extractFile(files(i));
%     [x{i},y{i},imageData{i}] = readOneTile(reg,num,pl,reso);
%     fprintf('tile %02d done\n',i)
% end
% save(['G:\BIGDATA\TOPIC 2\',sprintf('BirmLocalMap_reso%dm.mat',reso)],'x','y','imageData');



