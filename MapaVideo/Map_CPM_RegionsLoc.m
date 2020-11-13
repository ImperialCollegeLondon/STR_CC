figure;
hold on;

UKMap = getUKMap();

set(gca,'Linewidth',2);
Z = cpm*sf;
hh = pcolor(E,N,Z);shading flat
% cptcmap('precip_annualMeanUK', 'mapping','direct');%,'ncol',20);caxis([0,3100]);
cptcmap('BrBG_11', 'mapping','direct','ncol',9)
% axis off
caxis([300,2100])%([1,5.5]*365)
c = colorbar('location','Manual', 'position', [0.65 0.6 0.04 0.3],'fontsize',10);

c.Ruler.TickLabelFormat='%gmm';
c.FontSize = 8;%10;
c.Ruler.TickLength = [0.04,0.2];
ha = plot(UKMap.borderE/1000,UKMap.borderN/1000,'-','linewidth',1,'color',[0.5,0.5,0.5]);
ylim([0,1200]);xlim([0,800]);

% scalebar(gca,'Unit','km')

cmap = cptcmap('flood_blue','mapping','scaled','ncol',5,'flip',true);

angle = 15;

regionName = 'SUK';
leng = 88;
ul_x = 50;ul_y = 250;
ur_x = 660;ur_y = ul_y+tan(angle/180*pi)*(ur_x-ul_x);
X = [ul_x,ul_x+leng,ur_x+leng,ur_x,ul_x];
Y = [ul_y,ul_y-(leng/tan(angle/180*pi)),ur_y-(leng/tan(angle/180*pi)),ur_y,ul_y];
oneRegion(X,Y,regionName,cmap(1,:))
%
regionName = 'NWUK';
ll_x = X(1);ll_y = Y(1);
lr_x = 400;lr_y = ll_y+tan(angle/180*pi)*(lr_x-ll_x);
leng = 200;
X = [ll_x-leng,ll_x,lr_x,lr_x-leng,ll_x-leng];
Y = [ll_y+(leng/tan(angle/180*pi)),ll_y,lr_y,lr_y+(leng/tan(angle/180*pi)),ll_y+(leng/tan(angle/180*pi))];
oneRegion(X,Y,regionName,cmap(2,:))
%
regionName = 'NEUK';
ll_x = X(3);ll_y = Y(3);
lr_x = 600;lr_y = ll_y+tan(angle/180*pi)*(lr_x-ll_x);
leng = 200;
X = [ll_x-leng,ll_x,lr_x,lr_x-leng,ll_x-leng];
Y = [ll_y+(leng/tan(angle/180*pi)),ll_y,lr_y,lr_y+(leng/tan(angle/180*pi)),ll_y+(leng/tan(angle/180*pi))];
oneRegion(X,Y,regionName,cmap(3,:))
%

hold off
axis equal
xlim([0,800]);

ax = gca;
ax.TickDir = 'out';
ax.XTick = [];
ax.YTick = [];
ax.LineWidth = 0.5;
box on

hold on;
xloc = 100;
yloc = 1150;
barLength = 200;
tickLength = 30;
smalltickLength = 20;
setScaleBar(xloc,yloc,barLength,tickLength,smalltickLength);
%
savePath = 'E:\OneDrive - Imperial College London\dropbox\Fig_UKCP\';
fileName = sprintf('StudyRegion');
savePlot([savePath,filesep,fileName],'targetSize','1c','needreply','Y','onlyPng',false);

%%
function setScaleBar(xloc,yloc,barLength,tickLength,smalltickLength);
plot([xloc; xloc+barLength], ones(2,1)*yloc, '-k','LineWidth',1);hold on
plot([0,1/2,1].*[1,1,NaN]'*barLength+xloc,...
    [yloc-tickLength/2;yloc+tickLength/2;NaN],'k-','LineWidth',0.5);
plot(linspace(0,1,10-1).*[1,1,NaN]'*barLength+xloc,...
    [yloc-smalltickLength/2;yloc+smalltickLength/2;NaN],'k-','LineWidth',0.5);
ylim([0,1200]);xlim([0,800]);
fontsize = 7;
text(xloc,yloc, '0', 'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','fontsize',fontsize)
text(xloc+barLength/4,yloc, '50', 'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','fontsize',fontsize)
text(xloc+barLength/2,yloc, '100', 'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','fontsize',fontsize)
text(xloc+barLength,yloc, '200', 'HorizontalAlignment','center',...
    'VerticalAlignment','bottom','fontsize',fontsize)
text(xloc+barLength+30,yloc, 'km', 'HorizontalAlignment','left',...
    'VerticalAlignment','bottom','fontsize',fontsize)
end


function oneRegion(X,Y,regionName,color)
line(X,Y,'color',color);
hold on;
[x0,y0] = centroid(polyshape({X},{Y}));
text(x0,y0,regionName,...
    'fontsize',12,'horizontalalignment','center','fontweight','bold','background',[0.9 0.9 0.9 0.6]);
end


