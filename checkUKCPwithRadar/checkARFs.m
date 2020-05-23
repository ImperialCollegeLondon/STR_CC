% This file is to compute the Areal Reduction Factors (ARFs) for several datasets
% Method: Fixed Area Method
% Reference: (Svensson, *2010*)
%
% @ Yuting Chen
% Imperial College London
% yuting.chen17@imperial.ac.uk


%% CONFIGURATION

clear;clc;
close all
startUp();

REGIONS = REGIONS_info();
region = REGIONS.UK;


%% GET DATA TIME CONSUMING
%{
datasource = 'rad';
RES = getARFs_TS(datasource);
save('D:\UKCP18\ARFs\rad.mat','RES','datasource');

datasource = 'gear_1hr';
RES = getARFs_TS(datasource);
save('D:\UKCP18\ARFs\gear_1hr.mat','RES','datasource');


datasource = 'cpm2.2';
RES = getARFs_TS(datasource);
save('D:\UKCP18\ARFs\cpm.mat','RES','datasource');
%}

%% COMPUTE ARFs
figure;
setFigureProperty('Paper');
h = [];

DH = [1,2,6,24,48,25*24];
DH = [1,6,24];

simCol = [1 0 0];
obsCol = [0.3 0.3 0.3];
hind = 1;

ha = tight_subplot(1,3,[.1 .05],[.2 .15],[.1 .05]);
set(gcf,'units','points','position',[150,0,700,500]);


for dh = DH
    
    axes(ha(hind))
    
    datasource = 'cpm2.2';
    filename = 'D:\UKCP18\ARFs\cpm.mat';
    ARFs_sim = computeARFs(datasource,filename,dh);
    
    datasource = 'rad';
    filename = 'D:\UKCP18\ARFs\rad.mat';
    ARFs_rad = computeARFs(datasource,filename,dh);
    
    datasource = 'gear_1hr';
    filename = 'D:\UKCP18\ARFs\gear_1hr.mat';
    ARFs_gear = computeARFs(datasource,filename,dh);
    
    % PLOT RESULT
    %%
    [val,xarea] = deal([]);
    for ensNo = 1:12
        arfs_sim_thisEns = ARFs_sim(ARFs_sim.ensNo==ensNo,:);
        [val0,xarea0] = deal([]);
        for pind = 1:size(arfs_sim_thisEns,1)
            val0 = [val0;[1,nanmedian(arfs_sim_thisEns.arf{pind,1},1)]];
            xarea0 = [xarea0;[2.2.^2,arfs_sim_thisEns.trueArea2_2{pind}]];
            % plot(xarea(end,:),val(end,:),'.:','color',cmap(hind,:),'linewidth',1);hold on
        end
        val = cat(3,val,reshape(val0,[size(val0),1]));
        xarea = cat(3,xarea,reshape(xarea0,[size(xarea0),1]));
    end
    xarea = squeeze(nanmean(nanmean(xarea,1),3));
    val = squeeze(nanmean(val,1));
    hsim = plot(xarea,nanmedian(val,2),'-','color',[simCol,0.5],...
        'markerfacecolor',[simCol],'markersize',3);hold on;%'o-'
    sm = @(x)reshape(smooth(x,1),[1,numel(x)]);
    hsimran = fill([xarea,flip(xarea)],[sm(nanmin(val,[],2)),sm(flip(nanmax(val,[],2)))],...
        simCol,'LineStyle','none');
    alpha(0.3);
    
    [val,xarea] = deal([]);
    for pind = 1:size(ARFs_rad,1)
        val = [val;[1,nanmedian(ARFs_rad.arf{pind,1},1)]];
        xarea = [xarea;[2.2.^2,ARFs_rad.trueArea2_2{pind}]];
        % plot(xarea(end,:),val(end,:),'.--','color',cmap(hind,:),'linewidth',1);hold on
    end
    hrad = plot(nanmean(xarea,1),nanmean(val,1),'^--','color',[obsCol,0.5],...
        'linewidth',2,'markersize',3);
    hold on;
    
    [val,xarea] = deal([]);
    for pind = 1:size(ARFs_gear,1)
        val = [val;[1,nanmean(ARFs_gear.arf{pind,1},1)]];
        xarea = [xarea;[2.2.^2,ARFs_gear.trueArea2_2{pind}]];
        % plot(xarea(end,:),val(end,:),'.--','color',cmap(hind,:),'linewidth',1);hold on
    end
    hgear = plot(nanmean(xarea,1),nanmean(val,1),'o-','color',[obsCol,0.5],...
        'markerfacecolor',[obsCol]+0.5,'markersize',3);
    hold on;
    
    hind = hind+1;
    
    legend([hsim,hrad,hgear],{'cpm2.2','rad','gear'},'Location','southwest')
    title(strcat(string(dh),'h'))
    
    formatIt(dh)
    
end
% legend(h,{'60min','2h','6h','24h','2d','25d'})
format_xylabel(ha,1,3)

filePath = 'C:\Users\Yuting Chen\Dropbox (Personal)\Data_PP\Fig_UKCP\TemporalAnalysis\';
filename = [filePath,filesep,'ARFs_3durations'];
savePlot(filename,'units','centimeters','XYWH',[5,0,19,7],'needreply','N');


function formatIt(dh)
logsf=1;
xlim([2.2^2,1e4])
if dh<=24
    ylim([0.4,1])
else
    ylim([0.8,1])
end
grid minor
set(gca,'XScale','log')
xticks([1,10,100,1000,10000]*logsf)
xticklabels({'1','10','100','1000','10000'})
yticks([0.4,0.6,0.8,1]*logsf)
yticklabels({'0.4','0.6','0.8','1.0'})
set(gca,'linewidth',1)
% text(0,0.999,getSeasonName(thisseason),'fontsize',12);
% ax = gca;
% ax.FontSize = 20;
% ax.XTickLabelRotation = 90;
xlabel('Area (sq-km)');
ylabel('Areal Reduction Factors')
end

%% AXULLARY FUNCTION
function ARFs = computeARFs(datasource0,filename,dh)

load(filename,'RES','datasource');
% check correct datasource

if ~strcmpi(datasource,datasource0)
    error('Check datasource0 and filename');
end

try
    ARFs = table();
    for enNo = 1:size(RES{1},1)
        
        for pind = 1:numel(RES)
            res0 = RES{pind};
            [AM_p,AM_pa,arf,trueArea2_2] = deal([]);
            for yri = 1:size(res0,2)
                aggPoint = aggregate(res0(enNo,yri).PointTS,dh,'mean');
                pm_point = nanmax(aggPoint);
                pm_area = NaN(1,size(res0(enNo,yri).AreaTS,1));
                for ar = 1:size(res0(enNo,yri).AreaTS,1)
                    aggArea = aggregate(res0(enNo,yri).AreaTS(ar,:),dh,'mean');
                    pm_area(ar) = aggPoint(find(aggArea == nanmax(aggArea),1));
                end
                AM_p = [AM_p;pm_point];
                AM_pa = [AM_pa;pm_area];
                arf = [arf;pm_area./pm_point];
                
            end
            for ar = res0(enNo,yri).areaPerLen2_2
                ci = size(res0(enNo,1).isRainMat,1)/2;
                inds = ci-(ar-1)/2:ci+(ar-1)/2;
                trueArea2_2 = [trueArea2_2,...
                    nansum(reshape(res0(enNo,yri).isRainMat(inds,inds),[],1))*(2.2^2)];
            end
            area = (res0(enNo,yri).areaPerLen2_2*2.2).^2;
            if nansum(area-trueArea2_2)>1e-6
                fprintf('not error, but need check to confirm\n');
            end
            
            tab1 = table({AM_p},{AM_pa},{arf},{area},{trueArea2_2},enNo);
            ARFs = [ARFs;tab1];
        end
        
        
    end
    ARFs.Properties.VariableNames = {'AM_p','AM_pa','arf','area','trueArea2_2','ensNo'};
catch me
    me
end
end

function RES = getARFs_TS(datasource,options)

arguments
    datasource (1,:) char
    options (1,1) struct = struct('pNo',4,'period','historical');
end

RES = [];
[LOCS] = getLoc(options.pNo); % # <double>
% By now, LOCS are for SOUTH of EANGLAND.


areaPerLen = [1,3,5,11,15,21,45];

for pos = 1:size(LOCS,1)
    try
        region = struct('Name','temp','minE',LOCS(pos,1)-55,'minN',LOCS(pos,2)-55,'dimE',50,'dimN',50,'dx',2.2);
        isRain = 0;
        
        switch(datasource)
            
            case 'cpm2.2'
                YEARVEC = 1981:2000;
                ENSEMBLENO = getEnsNos();
                filePath = ['K:/UkCp18/',ENSEMBLENO{1},'/pr_rcp85_land-cpm_uk_2.2km_'];
                ncFileName = [filePath,ENSEMBLENO{1},'_1hr_19910801-19910830.nc'];
                LAT=ncread(ncFileName,'latitude');LON=ncread(ncFileName,'longitude');
                [E, N] = ll2os(LAT, LON); E = E/1000; N = N/1000;
            case 'rad'
                YEARVEC = 2007:2018;
                filePath = ['K:/UK_Radar_NetCDF/'];
                ncFileName = [filePath,'pr_nimrod_uk_2.2km_1hr_201809.nc'];
                E = ncread(ncFileName,'E'); N = ncread(ncFileName,'N');
                [N,E] = meshgrid(N,E);
                ENSEMBLENO = 1;
            case 'gear'
                
            case 'gear_1hr'
                YEARVEC = 1990:2014;
                filePath = 'K:/GEAR-1hr/';
                ncFileName = [filePath,'CEH-GEAR-1hr_199501.nc'];
                E = ncread(ncFileName,'x'); N = ncread(ncFileName,'y');%max to min
                [N,E] = meshgrid(N/1000,E/1000);
                ENSEMBLENO = 1;
            otherwise
                
                
        end
        [region.i,region.j] = getRegionIJ(E,N,region.minE,region.minN);
        warning off
        
        [AREA,POINT] = deal([]);
        
        for M=1:length(ENSEMBLENO)
            for year = YEARVEC
                
                [areal,point] = deal([]);
                tic
                for mon = 1:12
                    
                    if ~((strcmpi(datasource,'cpm2.2') && year==1980 && mon ~=12) ...
                            || (strcmpi(datasource,'cpm2.2') && year == 2000 && mon == 12))
                        
                        if strcmpi(datasource,'cpm2.2')
                            [filePath,fileName] = getFileName(datasource,year,mon,ENSEMBLENO{M});
                            listaRain = [filePath,fileName]; % A = ncinfo(listaRain);
                            try
                                rain_temp = squeeze(ncread(listaRain,'pr',...
                                    [region.i,region.j,1,1],[region.dimE,region.dimN,...
                                    Inf,...
                                    1]));
                            catch me
                                fprintf('Check:%s',listaRain);
                                me
                                rain_temp = NaN(region.dimE,region.dimN,24*30);
                            end
                        elseif strcmpi(datasource,'rad')
                            [filePath,fileName] = getFileName(datasource,year,mon);
                            listaRain = [filePath,fileName]; % A = ncinfo(listaRain);
                            rain_temp = squeeze(ncread(listaRain,'pr',...
                                [region.i,region.j,1],[region.dimE,region.dimN,Inf]));
                        elseif strcmpi(datasource,'gear_1hr')
                            % note: save format of gear is different from other
                            % dataset, where
                            % 1. northing vec is from max to min.
                            % 2. 1km resolution
                            [filePath,fileName] = getFileName(datasource,year,mon);
                            listaRain = [filePath,fileName]; %A = ncinfo(listaRain);
                            rain_temp = squeeze(ncread(listaRain,'rainfall_amount',...
                                [region.i,region.j+(-region.dimN)*region.dx+1,1],...
                                [region.dimE*region.dx,region.dimN*region.dx,Inf]));
                            rain_temp = rain_temp(:,end:-1:1,:);
                            % aggregate to 2.2
                            rain_temp = imresize3(rain_temp,'Scale',[5,5,1],'Method','nearest');
                            rain_temp = imresize3(rain_temp,'Scale',[1/11,1/11,1],'Method','box');
                        end
                        
                        
                        areal_temp = computeArea(rain_temp, areaPerLen);
                        point_temp = computepoint(rain_temp);
                        indRain = 1+size(areal,2):length(point_temp)+size(areal,2);
                        areal(:,indRain) = areal_temp;
                        point(indRain) = point_temp;
                        isRain = logical(isRain+logical(squeeze(nansum(rain_temp,3))));
                        
                    end
                    
                end
                toc
                fprintf('---%s-Year%02d---\n',region.Name,year);
                AREA{M,year-YEARVEC(1)+1} = areal;
                POINT{M,year-YEARVEC(1)+1} = point;
            end
        end
        switch(datasource)
            case 'cpm2.2'
                ARFs = struct('center',LOCS(pos,:),'AreaTS',AREA,'PointTS',POINT,...
                    'isRainMat',isRain,'areaPerLen2_2',areaPerLen,...
                    'year',YEARVEC);% 'Ensembles',ENSEMBLENO,
            case 'rad'
                ARFs = struct('center',LOCS(pos,:),'AreaTS',AREA,'PointTS',POINT,...
                    'isRainMat',isRain,'Product','Radar2_2','areaPerLen2_2',areaPerLen,...
                    'year',YEARVEC);
            case 'gear'
                
            case 'gear_1hr'
                ARFs = struct('center',LOCS(pos,:),'AreaTS',AREA,'PointTS',POINT,...
                    'isRainMat',isRain,'Product','GEAR_1hr','areaPerLen2_2',areaPerLen,...
                    'year',YEARVEC);
            otherwise
                
                
        end
        
        RES{pos} = ARFs;
        
    catch me
        me
        save('D:\UKCP18\ARFs\temp.mat','RES','LOCS','datasource','-v7.3')
    end
    
end

end

function [filePath,fileName] = getFileName(datasource,year,mon,No)
arguments
    datasource (1,:) char
    year (1,1) double
    mon (1,1) double
    No (1,:) char = ''
end

switch(datasource)
    case 'cpm2.2'
        filePath = ['K:/UkCp18/',No,'/'];
        fileName = sprintf('pr_rcp85_land-cpm_uk_2.2km_%s_1hr_%04d%02d01-%04d%02d30.nc',...
            No,year,mon,year,mon);
    case 'rad'
        filePath = 'K:/UK_Radar_NetCDF/';
        fileName = ['pr_nimrod_uk_2.2km_1hr_',sprintf('%04d%02d.nc',year,mon)];
    case 'gear_1hr'
        filePath = 'K:/GEAR-1hr/';
        fileName = ['CEH-GEAR-1hr_',sprintf('%04d%02d.nc',year,mon)];
end
end

function areal_temp = computeArea(rain_temp, areaPerLen);
i0 = round(size(rain_temp,1)/2);
areaPerLen = reshape(areaPerLen,[],1);
is = [i0-(areaPerLen-1)/2,i0+(areaPerLen-1)/2];
areal_temp = zeros(size(is,1),size(rain_temp,3));
for i = 1:size(is,1)
    thisArea = rain_temp(is(i,1):is(i,2),is(i,1):is(i,2),:);
    areal_temp(i,:) = squeeze(nansum(nansum((thisArea),1),2))./1;
end
end

function point_temp = computepoint(rain_temp);
i = round(size(rain_temp,1)/2);
point_temp = squeeze(rain_temp(i,i,:));
end


function [LOCS] = getLoc(pNo); % # <double>
% GETLOC() will pick up several points as well as corresponding rectangular area,
% with the larged area of 1000 0 sq-km

getPos = @(region)round([region.minE+50,region.minN+50]);

REGIONS = REGIONS_info();
LOCS = [];
LOCS(1,:) = getPos(REGIONS.EAng);
LOCS(2,:) = getPos(REGIONS.Westuk);
LOCS(3,:) = getPos(REGIONS.Birm);
LOCS(4,:) = getPos(REGIONS.London);
LOCS(5,:) = getPos(REGIONS.SWestuk);

end




