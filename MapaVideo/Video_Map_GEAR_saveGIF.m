clear;
clc;


%% Configuration
fold1='K:\GEAR-1hr\CEH-GEAR-1hr_199407.nc';
lista=dir(fold1);
StationID=fullfile({lista.folder},{lista.name});

ENSEMBLENO={'01','04','05','06','07','08','09','10','11','12','13','15'};

sizeRow = 50;
sizeCol = 50;
IntFac = 32; % to make the resolution = (Nimrod)

%% Let it Run

ncFileName = 'K:\GEAR-1hr\CEH-GEAR-1hr_199407.nc';
LAT=ncread(ncFileName,'lat');
LON=ncread(ncFileName,'lon');
[E, N] = ll2os(LAT, LON);
E = E/1000;
N = N/1000;

UKMap = load('H:\CODE_MATLAB\SpatialTemporalDATA\shapeFileFolder\UKBorderGrid.mat');


warning off
ntag = 1;


h = figure;
setFigureProperty('Paper');
XYWH = [50,50,400,600];
set(gcf,'units','points','position',XYWH,'defaultTextFontSize',20);
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'AnimatedGEARRain.gif';


listaRain='K:\GEAR-1hr\CEH-GEAR-1hr_201301.nc';

A=ncinfo(listaRain);

for time1 = 1:100%A.Variables(6).Dimensions(:,3).Length
    rr = ncread(listaRain,'rainfall_amount',...
        [1,1,time1],...
        [Inf,Inf,1]);
    % rr(rr<=1/32) = NaN;
    rr(rr==0) = NaN;
    pcolor(E,N,rr);
    shading flat;
    cb = colorbar;
    cptcmap('precip_meteoswiss', 'mapping','direct');%,'ncol',20);
    hold on;
    plot(UKMap.borderE/1000,UKMap.borderN/1000,'k-','linewidth',1);
    set(gca,'Visible','off');
    caxis([0,40])
    timn=datetime(2013,1,1)+time1/24;
    text(400,-10,sprintf('%s',timn),'color','k','fontweight', 'bold' )
    set(gca,'Color','none')
    drawnow
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if ntag == 1
        imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',0.2);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.2);
    end
    clf(h)
    ntag = ntag+1;
    if ntag>24*7
        continue;
    end
    
    
    
end




