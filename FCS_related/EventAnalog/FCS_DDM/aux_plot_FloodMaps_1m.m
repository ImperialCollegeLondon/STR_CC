%
% aux_plot_FloodMaps_1m()
%

%% load data
[floodmaps,E,N] = aux_temp_loadFloodMaps_WCs(90);
[fli,flj,S] = find(floodmaps{1}); 
S = 100*S;
%% determine central location to show

% scatter(E(i),N(j),S*50,'r','fill');
[pE,pN] = getLoc(E,N,fli,flj);


%% find out corresponding tile from OS local map

load('H:\CODE_MATLAB\SpatialTemporalDATA\shapeFileFolder\BirmLocalMap_reso2m.mat','x','y','imageData');
Area = getBoundaryShp('Birm');

f = figure;
f.GraphicsSmoothing = 'off';
ax1 = axes;
XLIM = [Inf,0];YLIM = [Inf,0];
for i = 1:length(x)
    for tileNum = 1:numel(x{i})
        
        if abs(x{i}{tileNum}-pE)<5000 & abs(y{i}{tileNum}-pN)<5000
            pcolor(x{i}{tileNum},y{i}{tileNum},imageData{i}{tileNum}); shading flat; hold on;
            XLIM = [min(XLIM(1),min(x{i}{tileNum})),max(XLIM(2),max(x{i}{tileNum}))];
            YLIM = [min(YLIM(1),min(y{i}{tileNum})),max(YLIM(2),max(y{i}{tileNum}))];
        else
            fprintf('Not this tile\n')
        end
        
    end
    cptcmap('GMT_gray','mapping', 'scaled','flip',false,'ncol',256); caxis([0,255])
    drawnow;
end
% plot(Area.X,Area.Y,'k','linewidth',2)
alpha(0.5);
set(gca,'YDir','normal')
axis equal
axis off
xlim(XLIM);ylim(YLIM)
clear x y imageData i j
toc

%% imshow background


%% scatter floodExtent
view(2)
ax2 = axes;

axislim = round(nanmax(S(:))/10);
[cmap, lims, ticks, bfncol, ctable] = cptcmap('flood_blue','mapping','scaled','ncol',15);
colfmap = NaN(numel(S(:)),3);
colfmap(~isnan(S(:)),:) = cmap(getLevel(S(~isnan(S)),linspace(0,axislim,15)),:);
scatter(ax2,E(fli),N(flj),8,colfmap,'fill');shading flat;


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

xlim(XLIM);ylim(YLIM)

%% AUXILLARY FUNC
function [pE,pN]  = getLoc(E,N,i,j)
pE = [414500];
pN = [292100];
end
function mapLevel = getLevel(maps,lowThre)
mapLevel = NaN(size(maps));
lowThre = [lowThre,inf];
for li = 1:numel(lowThre)-1
    mapLevel(maps>=lowThre(li) & maps<lowThre(li+1)) = li;
end
end
