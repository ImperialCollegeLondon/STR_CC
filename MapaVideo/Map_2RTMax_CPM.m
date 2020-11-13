figure;
hold on;
RT = 5;
load(['J:\D\UKCP18\UK\AM_CPM_1980-2000_RT',num2str(RT),'.mat'])
CPM_saved = load('J:\D\UKCP18\UK\MonMean_CPM.mat','E','N');
UKMap = getUKMap();

I = squeeze(nanmean(cell2mat(cellfun(@(x)reshape(x,[1,size(x)]),I','UniformOutput',false)),1));
in = getTrimTag('unit','km','product','cpm2.2');
I = conv2(I,ones(6,6)/36,'same');
I(~in) = NaN;

E = CPM_saved.E;N = CPM_saved.N;

pcolor(E,N,I);shading flat;hold on
cptcmap('GnBu_09', 'mapping','scaled','ncol',9);
caxis([20,30])
% caxis([prctile(I(:),5),prctile(I(:),99)])
axis('equal');
axis off
c = colorbar('location','Manual', 'position', [0.65 0.6 0.04 0.3],'fontsize',10);

c.Ruler.TickLabelFormat='%gmm/h';c.FontSize = 8;
plot(UKMap.borderE/1000,UKMap.borderN/1000,'-','linewidth',1,'color',[0.5,0.5,0.5]);
ylim([0,1200]);xlim([0,800]);
nanmean(I(:))

%
cmap = cptcmap('flood_blue','mapping','scaled','ncol',5,'flip',true);
angle = 15;

regionName = 'SUK';
leng = 88;
ul_x = 50;ul_y = 250;
ur_x = 660;ur_y = ul_y+tan(angle/180*pi)*(ur_x-ul_x);
X = [ul_x,ul_x+leng,ur_x+leng,ur_x,ul_x];
Y = [ul_y,ul_y-(leng/tan(angle/180*pi)),ur_y-(leng/tan(angle/180*pi)),ur_y,ul_y];
oneRegion(X,Y,regionName,cmap(1,:))

regionName = 'NWUK';
ll_x = X(1);ll_y = Y(1);
lr_x = 400;lr_y = ll_y+tan(angle/180*pi)*(lr_x-ll_x);
leng = 200;
X = [ll_x-leng,ll_x,lr_x,lr_x-leng,ll_x-leng];
Y = [ll_y+(leng/tan(angle/180*pi)),ll_y,lr_y,lr_y+(leng/tan(angle/180*pi)),ll_y+(leng/tan(angle/180*pi))];
oneRegion(X,Y,regionName,cmap(2,:))

regionName = 'NEUK';
ll_x = X(3);ll_y = Y(3);
lr_x = 600;lr_y = ll_y+tan(angle/180*pi)*(lr_x-ll_x);
leng = 200;
X = [ll_x-leng,ll_x,lr_x,lr_x-leng,ll_x-leng];
Y = [ll_y+(leng/tan(angle/180*pi)),ll_y,lr_y,lr_y+(leng/tan(angle/180*pi)),ll_y+(leng/tan(angle/180*pi))];
oneRegion(X,Y,regionName,cmap(3,:))

hold off
axis equal
ylim([0,1200]);xlim([0,800]);

ax = gca;
ax.XTick = [];
ax.YTick = [];
ax.LineWidth = 0.5;
box on
axis on

savePath =  'E:\OneDrive - Imperial College London\dropbox\Fig_UKCP\';
fileName = sprintf('AM%dyear',RT);
savePlot([savePath,filesep,fileName],'targetSize','1c','needreply','Y','onlyPng',false);

function oneRegion(X,Y,regionName,color)
line(X,Y,'color',color);
hold on;
[x0,y0] = centroid(polyshape({X},{Y}));
% text(x0,y0,regionName,'fontsize',12,'horizontalalignment','center','fontweight','bold','background',[0.9 0.9 0.9 0.6]);
end
