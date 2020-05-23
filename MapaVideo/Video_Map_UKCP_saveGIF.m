clear;
clc;


%% Configuration
fold1='C:/Data/UK/Midas/Stations/Final Files/*.nc';
lista=dir(fold1);
StationID=fullfile({lista.folder},{lista.name});

ENSEMBLENO={'01','04','05','06','07','08','09','10','11','12','13','15'};

sizeRow = 50;
sizeCol = 50;
IntFac = 32; % to make the resolution = (Nimrod)

REGIONS = REGIONS_info();
region = REGIONS.SWestuk;
% region = struct('Name','SWestuk','minE',221,'minN',41,'dimE',50,'dimN',50);
% region = struct('Name','Westuk','minE',219.8,'minN',225.8,'dimE',50,'dimN',50);
% region = struct('Name','Scotland','minE',206.7,'minN',627.6,'dimE',50,'dimN',50);
% region = struct('Name','London','minE',475,'minN',125.3,'dimE',50,'dimN',50);

%% Let it Run 
filePath = ['J:/UkCp18/',ENSEMBLENO{1},'/pr_rcp85_land-cpm_uk_2.2km_'];
ncFileName = [filePath,ENSEMBLENO{1},'_1hr_19910801-19910830.nc'];
LAT=ncread(ncFileName,'latitude');
LON=ncread(ncFileName,'longitude');
[E, N] = ll2os(LAT, LON);
E = E/1000;
N = N/1000;

UKMap = load('H:\CODE_MATLAB\SpatialTemporalDATA\Radar\UKBorderGrid.mat');

warning off
ntag = 1;


h = figure;
XYWH = [50,50,500,500];
set(gcf,'units','points','position',XYWH);
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'AnimatedUKCPRain.gif';

for mon = 1:12
    
    RainEnsembles = cell(length(ENSEMBLENO),1);
    
    mkdir(['F:/UKCP18/',region.Name]);
    saveDir = ['F:/UKCP18/',region.Name,'/Ensems',...
        '_mon',sprintf('%02d',mon)];
    
    for M=1:length(ENSEMBLENO)
        
        filePath = ['J:/UkCp18/',ENSEMBLENO{M},'/pr_rcp85_land-cpm_uk_2.2km_'];
        
        listaRain=dir([filePath,ENSEMBLENO{M},'_1hr_*',sprintf('%02d',mon),'30.nc']);
        
        listaRain=fullfile({listaRain.folder},{listaRain.name});
        
        Rain = int16(NaN(region.dimE,region.dimN,length(listaRain)*30*24)*length(ENSEMBLENO));
        
        for L=1:length(listaRain)
            
            A=ncinfo(listaRain{L});
            
            for time1 = 1:A.Variables(1).Dimensions(:,3).Length
                rr = ncread(listaRain{L},'pr',...
                    [1,1,time1,1],...
                    [Inf,Inf,1,...
                    A.Variables(1).Dimensions(:,4).Length]);
                % rr(rr<=1/32) = NaN;
                pcolor(E,N,rr);
                shading flat;
                cb = colorbar;
                cptcmap('precip_meteoswiss', 'mapping','direct');%,'ncol',20);
                hold on;
                plot(UKMap.borderE/1000,UKMap.borderN/1000,'w-','linewidth',1);
                set(gca,'Visible','off');
                caxis([0,40])
                timn=datetime(1981,1,1)+time1/24;
                text(490,-200,sprintf('%s',timn),'color','w','fontweight', 'bold' )
                set(gca,'Color','none')
                drawnow
                frame = getframe(h);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                if ntag == 1
                    imwrite(imind,cm,filename,'gif','Loopcount',inf,'DelayTime',0.1);
                else
                    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
                end
                clf(h)
                ntag = ntag+1;
                if ntag>24*7
                    continue;
                end
                
            end
            
        end
        
    end
    
    fprintf('---%s-Mon%02d---\n',region.Name,mon);
end
