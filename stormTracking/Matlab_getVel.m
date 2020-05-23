ENSEMBLENO = getEnsNos();

REGIONS = REGIONS_info();
IntFac = 32; % to make the resolution = (Nimrod)
options = 'Save';
for ENSEMBLENO_1 = ENSEMBLENO
    for MON = 6:8
        for region = [REGIONS.SCO,REGIONS.WAL,REGIONS.EUK]%[REGIONS.SWestuk,REGIONS.Westuk,REGIONS.Scotland,REGIONS.London]%
            
            data = struct('Years',[1980,2000],...
                'fileGetPath','K:/UkCp18/',...
                'savePath',['D:/UKCP18/',region.Name]);
            readCPM_pr(region,ENSEMBLENO_1,MON,IntFac,data,options);
            
            data = struct('Years',[2020,2040],...
                'fileGetPath','K:/UkCp18_FutureTemp/',...
                'savePath',['D:/UKCP18_Future/',region.Name]);
            readCPM_pr(region,ENSEMBLENO_1,MON,IntFac,data,options);
            
            data = struct('Years',[2060,2080],...
                'fileGetPath','K:/UkCp18_FutureTemp/',...
                'savePath',['D:/UKCP18_Future_2060_2080/',region.Name]);
            readCPM_pr(region,ENSEMBLENO_1,MON,IntFac,data,options);
            
        end
    end
end

%%
fprintf('EXTRACTION Done!\n')
UKCP_process_ConvStorm_JJA


%%
load('D:\UKCP18\bigSCO\Ensems_mon01.mat')
R = double(RainEnsembles{1}(:,:,200:300))/IntFac;
R = permute(R,[3,1,2]);
save('D:\UKCP18\bigSCO\oneM.mat','R')
%%

load('D:\UKCP18\bigSCO\oneM.mat','R')
zeroVal = 0;
R0 = R(10:30,:,:);
R = func_onlyMCS(R0,2,zeroVal);
R = func_R2dBZ(R,'UKMO',zeroVal);
mInt = nanmax(reshape(R0,[size(R0,1),prod(size(R0,[2,3]))]),[],2);

trajec = [];



for ts = 1:size(R0,1)
    Z = squeeze(R0(ts,:,:));
    s  = regionprops('table',logical(Z),mat2gray(Z,[0 65]),{'WeightedCentroid','Centroid','Area'});
    centroids = cat(1, s.WeightedCentroid);
    % centroids(s.Area<200)=[];
    trajec(ts,:) = [centroids(s.Area == max(s.Area),1), centroids(s.Area == max(s.Area),2)];
    Z = squeeze(R0(ts,:,:));
    Z(Z==zeroVal) = NaN;
    pcolor(Z);shading flat;alpha(0.2)
    cptcmap('precip_meteoswiss', 'mapping','direct');%,'ncol',20);
    hold on
    
    plotArgs = [];
    for i = 1 : size(trajec,1)-1
        thisInt = mean(mInt(i:i+1));
        colorVector = getCol(thisInt);
        line('XData', trajec(i:i+1,1), 'YData', trajec(i:i+1,2), ...
            'Color', colorVector,'linewidth',2);
    end
    
    % h = plot(trajec(:,1),trajec(:,2),'wo-','markerfacecolor','y','linewidth',2);
    
    ax = gca;
    ax.Color = 0.5*ones(1,3);
    colorbar
    hold off
    pause(0.5)
end

%%
zeroVal = 0;
load('D:\UKCP18\bigSCO\oneM.mat','R')
R = R(10:30,:,:);
R = func_onlyMCS(R,5,zeroVal);
R = func_R2dBZ(R,'UKMO',zeroVal);

R1 = squeeze(R(11,:,:));
R2 = squeeze(R(12,:,:));
save('Rainfall_dBZ_2','R1','R2');
%%
imagesc(R1)
cptcmap('precip_meteoswiss', 'mapping','scaled');%,'ncol',20);

figure
imagesc(R2)
cptcmap('precip_meteoswiss', 'mapping','scaled');%,'ncol',20);

%%

zeroVal = 0;
load('D:\UKCP18\bigSCO\oneM.mat','R')
R = R(10:30,:,:);

vidReader = VideoReader('visiontraffic.avi','CurrentTime',11);

opticFlow = opticalFlowLK('NoiseThreshold',0.009);
h = figure;
movegui(h);
hViewPanel = uipanel(h,'Position',[0 0 1 1],'Title','Plot of Optical Flow Vectors');
hPlot = axes(hViewPanel);

for sni = 1:size(R,1)
    
    frameGray = squeeze(R(sni,:,:));
    frameGray = func_onlyMCS(frameGray,5,zeroVal);
    frameGray = func_R2dBZ(frameGray,'UKMO',zeroVal);
    % frameGray = func_conv(frameGray,5);
    flow = estimateFlow(opticFlow,frameGray);
    vx = flow.Vx(frameGray ~= zeroVal);
    vy = flow.Vy(frameGray ~= zeroVal);
    
    frameGray(frameGray==zeroVal) = NaN;
    pcolor(frameGray);shading flat;
    cptcmap('precip_meteoswiss', 'mapping','scaled');%,'ncol',20);
    caxis([35,50])
    hold on
    hh = plot(flow,'DecimationFactor',[5 5],'ScaleFactor',10,'Parent',hPlot);
    hold off
    % pause(10^-1)
    axis normal
    
end


function colvec = getCol(val)
[cmap, lims, ticks, bfncol, ctable] = cptcmap('GMT_no_green','mapping','scaled','ncol',15);
colvec = cmap(getLevel(val,linspace(0,10,15)),:);
    function mapLevel = getLevel(val,lowThre)
        mapLevel = NaN(size(val));
        lowThre = [lowThre,inf];
        for li = 1:numel(lowThre)-1
            mapLevel(val>=lowThre(li) & val<lowThre(li+1)) = li;
        end
    end
end

function R = func_onlyMCS(R,thres,zeroVal);
R(R<thres) = zeroVal;
end

function R = func_R2dBZ(R,source,zeroVal)
if strcmp(source,'UKMO')
    R(R~=zeroVal) = (log10(200)+1.6*log10(R(R~=zeroVal)))*10;
end
end

function R = func_dBZ2R(R,source,zeroVal)
if strcmp(source,'UKMO')
    R(R~=zeroVal) = ((10.^(R(R~=zeroVal)/10))/200).^(5/8);
end
end

function R = func_conv(R,nn)
R = conv2(R,ones(nn,nn)/(nn*nn),'same');
end

