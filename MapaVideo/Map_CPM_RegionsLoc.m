figure;
hold on;

UKMap = getUKMap();

set(gca,'Linewidth',2);
Z = cpm*sf;
pcolor(E,N,Z);shading flat
cptcmap('precip_annualMeanUK', 'mapping','direct');%,'ncol',20);

% axis off

c = colorbar('location','Manual', 'position', [0.65 0.6 0.04 0.3],'fontsize',10);

c.Ruler.TickLabelFormat='%gmm';
c.FontSize = 10;
plot(UKMap.borderE/1000,UKMap.borderN/1000,'-','linewidth',1,'color',[0.5,0.5,0.5]);
ylim([0,1200]);caxis([0,3100]);xlim([0,800]);

cmap = cptcmap('flood_blue','mapping','scaled','ncol',5,'flip',true);

angle = 15;

regionName = 'SUK';
leng = 88;
ul_x = 50;ul_y = 250;
ur_x = 660;ur_y = ul_y+tan(angle/180*pi)*(ur_x-ul_x);
X = [ul_x,ul_x+leng,ur_x+leng,ur_x,ul_x];
Y = [ul_y,ul_y-(leng/tan(angle/180*pi)),ur_y-(leng/tan(angle/180*pi)),ur_y,ul_y];
oneRegion(X,Y,regionName,cmap(1,:))
%%
regionName = 'NWUK';
ll_x = X(1);ll_y = Y(1);
lr_x = 400;lr_y = ll_y+tan(angle/180*pi)*(lr_x-ll_x);
leng = 200;
X = [ll_x-leng,ll_x,lr_x,lr_x-leng,ll_x-leng];
Y = [ll_y+(leng/tan(angle/180*pi)),ll_y,lr_y,lr_y+(leng/tan(angle/180*pi)),ll_y+(leng/tan(angle/180*pi))];
oneRegion(X,Y,regionName,cmap(2,:))
%%
regionName = 'NEUK';
ll_x = X(3);ll_y = Y(3);
lr_x = 600;lr_y = ll_y+tan(angle/180*pi)*(lr_x-ll_x);
leng = 200;
X = [ll_x-leng,ll_x,lr_x,lr_x-leng,ll_x-leng];
Y = [ll_y+(leng/tan(angle/180*pi)),ll_y,lr_y,lr_y+(leng/tan(angle/180*pi)),ll_y+(leng/tan(angle/180*pi))];
oneRegion(X,Y,regionName,cmap(3,:))
%%
hold off
axis equal
xlim([0,800]);

ax = gca;
ax.XTick = [];
ax.YTick = [];
ax.LineWidth = 0.5;
box on

savePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\';
fileName = sprintf('StudyRegion');
savePlot([savePath,filesep,fileName],'units','centimeters','XYWH',[5,0,8,11],'needreply','Y','onlyPng',false);

%%
function oneRegion(X,Y,regionName,color)
line(X,Y,'color',color);
hold on;
[x0,y0] = centroid(polyshape({X},{Y}));
text(x0,y0,regionName,...
    'fontsize',12,'horizontalalignment','center','fontweight','bold','background',[0.9 0.9 0.9 0.6]);
end


