% CHECK SPATIAL CLUSTERING

clear;clc
close all

REGIONS = REGIONS_info();
region = REGIONS.smallUK;
ENSEMBLENO=getEnsNos();
RadVec = [1:4,[1:10]*5,60,75,100];

imageNo = 1;

for enNo = 1
    itag = 1;
    
    in = getTrimTag('unit','km','product','radar2.2');
    % in = getTrimTag('unit','km','product','cpm2.2');
    
    [AK_NC,AK] = deal([]);
    
    for year = 2007:2018
        Rain = [];
        
        tic
        for mon = 6:8
            try
                [E,N,Rain0,scaleF,region] = readRADAR_nc_1(region,year,mon,[],imageNo);
                ins = repmat(in,1,1,size(Rain0,3));
                Rain0(~ins) = NaN;
                region = REGIONS.SUK;
                [region.i,region.j] = getRegionIJ(E,N,region.minE,region.minN);
                Rain0 = Rain0(region.i:region.i+region.dimE-1,...
                    region.j:region.j+region.dimN-1,:);
                Rain = cat(3,Rain,Rain0);
                
            catch me
                fprintf('Issue!\n')
            end
        end
        toc
        
        pmax = squeeze(nanmax(nanmax(Rain,[],1),[],2));
        [~,I] = sort(pmax,'descend');
        intenseInd = I(1:24);
        % intenseInd = find(nansum(reshape(Rain,[],size(Rain,3))>50,1)>=2);
        % ind(nanmean(reshape(Rain,[],size(Rain,3)),1)>1);
        
        Rain(isnan(Rain)) = 0;
        Rain(Rain<0.1)=0;
        
        
        for ind = 1:numel(intenseInd)
            
            OneSnap = squeeze(Rain(:,:,intenseInd(ind)));
            OneSnap = imresize(OneSnap,'Scale',[0.5,0.5],'method','box');
            [K,K_NC,Q01,Q99] = deal([]);
            
            tic
            for radInd = 1:length(RadVec)
                [K(radInd),K_NC(radInd),Q01(radInd),...
                    Q99(radInd)] = RipleysK2(OneSnap,RadVec(radInd));
                fprintf('Radius %02d km Done\n',RadVec(radInd));
                % Q99(radInd)] = RipleysK1(Rain,RadVec(radInd));
            end
            toc
            
            plot(RadVec,sqrt(K/pi)-RadVec,'r-','linewidth',1);
            hold on;
            pause(0.02)
            AK_NC(itag,:) = K_NC;
            AK(itag,:) = K;
            itag = itag+1;
            
        end
        
        save(['D:/UKCP18/Radar_4d4_2007_',sprintf('%04d',year),'_JJA_pJJAp2_RipletK2.mat'],...
            'AK_NC','AK','RadVec');
    end
    grid minor
    xlabel('r(Km)');
    ylabel('K(r)');
    
    figure;
    plot(RadVec*2.2,median(AK_NC,1),'r-','linewidth',2);
    close all
    
    save(['D:/UKCP18/Radar_4d4_JJA_pJJAp2_RipletK2.mat'],'AK_NC','AK','RadVec');
end

% hold on;
% plot(RadVec,Q01,'r--',RadVec,Q99,'r--')

% read CPM one snapshot.
function [E,N,Rain,scaleF,region] = readRADAR_nc_1(region,year,mon,ensNo,imageNo)

fileName = sprintf('pr_nimrod_uk_2.2km_1hr_%04d%02d.nc',...
    year,mon);
filePath = ['K:/UK_Radar_NetCDF/'];
listaRain = [filePath,fileName];
A = ncinfo(listaRain);


ncFileName = [filePath,'pr_nimrod_uk_2.2km_1hr_201812.nc'];
E=ncread(ncFileName,'E');
N=ncread(ncFileName,'N');
[N,E] = meshgrid(N,E);
% E = E/1000;
% N = N/1000;

Rain = ncread(listaRain,'pr',...
    [1,1,imageNo],...
    [Inf,Inf,Inf]);
scaleF = 1;


end