% CHECK SPATIAL CLUSTERING

clear;clc
close all

REGIONS = REGIONS_info();
region = REGIONS.smallUK;
ENSEMBLENO=getEnsNos();
RadVec = [1:4,[1:10]*5,60,75,100];


imageNo = 1;
for enNo = 1:length(ENSEMBLENO)
    itag = 1;
    
    in = getTrimTag('unit','km','product','cpm2.2');
    ins = repmat(in,1,1,720);
    
    [AK_NC,AK] = deal([]);
    
    for year = 1981:2000
        Rain = [];
        for mon = 6:8
            try
                [E,N,Rain0,scaleF,region] = readCPM_nc_1(region,year,mon,ENSEMBLENO{enNo},imageNo);
                Rain0(~ins) = NaN;
                region = REGIONS.SUK;
                [region.i,region.j] = getRegionIJ(E,N,region.minE,region.minN);
                Rain0 = Rain0(region.i:region.i+region.dimE-1,...
                    region.j:region.j+region.dimN-1,:);
                
                % aggregate to 1day resolution
                Rain0 = convn(Rain0,ones(1,1,24),'valid');
                
                Rain = cat(3,Rain,Rain0);
            catch me
                fprintf('Issue!\n')
            end
        end
        pmax = squeeze(nanmax(nanmax(Rain,[],1),[],2));
        [~,I] = sort(pmax,'descend');
        I(I<I(1)+24 & I>I(1)-24 & I~=I(1)) = [];
        intenseInd = I(1:2); % I excluded those 'same' event/day, so now only can be 2%
        
        % intenseInd = find(nansum(reshape(Rain,[],size(Rain,3))>50,1)>=2);
        % ind(nanmean(reshape(Rain,[],size(Rain,3)),1)>1);
        
        Rain(isnan(Rain)) = 0;
        Rain(Rain<0.1)=0;
        
        for ind = 1:numel(intenseInd)
            
            OneSnap = squeeze(Rain(:,:,intenseInd(ind)));
            % aggregate to 4.4km for a quicker computing!
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
            
            plot(RadVec*4.4,sqrt(K/pi)-RadVec,'r-','linewidth',1);
            hold on;
            pause(0.02)
            AK_NC(itag,:) = K_NC;
            AK(itag,:) = K;
            itag = itag+1;
            
        end
        
    end
    grid minor
    xlabel('r(Km)');
    ylabel('K(r)');
    
    figure;
    plot(RadVec*4.4,median(AK_NC,1),'r-','linewidth',2);
    close all
    
    save(['D:/UKCP18/CPM_4d4_Ensno',ENSEMBLENO{enNo},'_JJA_pJJADaily_p2_RipletK2.mat'],...
        'AK_NC','AK','RadVec');
end


run('H:\CODE_MATLAB\SpatialTemporalDATA\checkUKCPwithRadar\checkCluster_Radar_Daily.m')

% hold on;
% plot(RadVec,Q01,'r--',RadVec,Q99,'r--')

% read CPM one snapshot.
function [E,N,Rain,scaleF,region] = readCPM_nc_1(region,year,mon,ensNo,imageNo)

fileName = sprintf('pr_rcp85_land-cpm_uk_2.2km_%s_1hr_%04d%02d01-%04d%02d30.nc',...
    ensNo,year,mon,year,mon);
filePath = ['K:/UkCp18/',ensNo,'/'];
listaRain = [filePath,fileName];
A = ncinfo(listaRain);


filePath = ['K:/UkCp18/01/pr_rcp85_land-cpm_uk_2.2km_'];
ncFileName = [filePath,'01','_1hr_19910801-19910830.nc'];
LAT=ncread(ncFileName,'latitude');
LON=ncread(ncFileName,'longitude');
[E, N] = ll2os(LAT, LON);
E = E/1000;
N = N/1000;

Rain = squeeze(ncread(listaRain,'pr',...
    [1,1,imageNo,1],...
    [Inf,Inf,Inf,1]));
scaleF = 1;


end