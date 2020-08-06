% LOAD DATA
expNo = 'C1P5In-Rad1';
% A = randi([1,157],1,5);
for eventNo4val = [48,61,90,94,103]
eviPlot = eventNo4val;
[version,testConfig,filefolder] = getExpNoInfo(expNo);
load([filefolder,sprintf('STATS_%s_ev%03d',version,eventNo4val)],'STATS','FOREOUTPUT');

E = FOREOUTPUT.E;
N = FOREOUTPUT.N;

floodmaps_test = squeeze(FOREOUTPUT.FLO_test.floodmaps(:,:,:));
floodmaps_pred = squeeze(FOREOUTPUT.FLO_pred.floodmaps(:,:,:,:));

% PLOT FLOODMAP
floodobs = floodmaps_test;
floodpre = floodmaps_pred;

savePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_FCS';
fileTag = sprintf('Web_evi%03d',eviPlot);
filename = [savePath,filesep,'Prediction_FloodMap_',fileTag,'_PredPrct50'];
Time = Times(eventNo4val);

save(filename,'floodobs','floodpre','E','N','Time')

floodpre = prctile(floodmaps_pred,50,4);
axislim = nanmax(floodobs(:));
ax1 = subplot(2,3,[1,2,3]);
plotFlodeoth(ax1,axislim,squeeze(prctile(floodobs,50,1)),E,N);

ax1 = subplot(2,3,[4]);
plotFlodeoth(ax1,axislim,squeeze(prctile(floodpre,10,1)),E,N);

ax1 = subplot(2,3,[5]);
plotFlodeoth(ax1,axislim,squeeze(prctile(floodpre,50,1)),E,N);

ax1 = subplot(2,3,[6]);
plotFlodeoth(ax1,axislim,squeeze(prctile(floodpre,90,1)),E,N);

pause(2)

end
%%



function [] = plotFlodeoth(ax1,axislim,floodval,E,N)

aggScale = 1; mode = 'max';
maps = aggregateImage(floodval,aggScale,mode);
E0 = imresize(E,'scale',1/aggScale,'method','box');
N0 = imresize(N,'scale',1/aggScale,'method','box');
plot_DemoFloodMaps(maps,E0,N0,axislim)

    function plot_DemoFloodMaps(maps,E0,N0,axislim)
        maps(maps == 0) = NaN;
        % pcolor(ax2,E,N,maps);shading flat;
        
        [cmap, lims, ticks, bfncol, ctable] = cptcmap('flood_blue','mapping','scaled','ncol',15);
        colfmap = NaN(numel(maps(:)),3);
        colfmap(~isnan(maps(:)),:) = cmap(getLevel(maps(~isnan(maps)),linspace(0,axislim,15)),:);
        scatter(E0(:),N0(:),(maps(:))*4,colfmap,'fill');shading flat;
        alpha(0.6)
       

        cptcmap('flood_blue', 'mapping','scaled','ncol',15);
        % axis off
        box on
        axis equal
        caxis([0,axislim]);
%         c = colorbar;
%         t0 = c.TickLabels;
%         c.TickLabels = strcat(t0,'cm');
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


