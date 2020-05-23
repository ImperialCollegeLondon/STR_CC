% -------------------------------------------------------%
% This file is to convert .mat day2day radar info. to netCDF4 file
% Data required:
% Radar Data: 
%            hourly radar info for whole uk save in '.mat'
%            (aggregated from 5min composite radar, using 'yt_func')
% Region Data:
%            can be any area within UK boundary
%            save in 'REGIONS_info.m'
% NOTICE: 
%            input NIMROD radar is <int16> format with a SF of 32.
% that means the resolution of aggregated(nansum) hourly radar is '1/32' mm/hour.
%
% Matlabe version:
% R2019b required.
%
% Output NC file info:
% pr: <single> [E,N,T]
% resolution: 2.2km, 1hour
% no scaleF
%
% @ Yuting Chen
% Update: 2019.01.07
% -------------------------------------------------------%


clear;clc
close all
addpath('C:\Users\Yuting Chen\Dropbox (Personal)\MatLab_Func');


%% CONFIGURATE INPUT
REGIONS_info = REGIONS_info();
RD_fp = 'H:\DATA_RADAR\UK_Radar_HourlyAggregate\';
RADAR_info = struct('row',2175,'col',1725);
[RADAR_info.E,RADAR_info.N,RADAR_info.EE,RADAR_info.NN] = ...
    getNimrodEN();% 'row',2175,'col',1725

Year1 = 2018;
Year2 = 2018;


%% Aggregate + Convert + Save Whole UK
%
% Size: 1GB per month
%
sfPath = 'K:\UK_Radar_NetCDF\';
mkdir(sfPath)
region = REGIONS_info.UK;
extractWholeUK(region,RADAR_info,RD_fp,Year1,Year2,sfPath)



function extractWholeUK(region,RADAR_info,RD_fp,Year1,Year2,sfPath)
%% LOAD and PROCESS data

TS = [];
scaleF = 32;
[loci,locj] = getLoc_IRADAR(RADAR_info,region);

for year = Year1:1:Year2
    for mon = 1:12
        
        PRS2_2 = NaN(region.dimN,region.dimE,0);% here, notice: NIMROD order is [N,E,T];
        for day = 1:eomday(year,mon)
		
		    RD_fn = [sprintf('%04d_%d_%d.mat',year,mon,day)];
			load([RD_fp,RD_fn]);
			[PRS0,RADAR_info] = getPRS(PRS0,RADAR_info,region,loci,locj);
            PRS0(PRS0 >= 100*scaleF | PRS0 < 0) = NaN;
            % threshold was chosen because maximum hourly rain:92mm in 12
            % July 1901 (MetOffice)
            % ref: https://www.metoffice.gov.uk/research/climate/maps-and-data/uk-climate-extremes
            
            % #imresize#
            % use Bilinear method, which is an interpolation which considers nearest 2-by-2 neighbourhood;
            % The method is choosen because a resolution of 2.2 km is close to the doule value of radar resolution (1km);
            % prs_temp = imresize3(PRS0,[region.dimN,region.dimE,size(PRS0,3)],'linear');
            
            % #imresize3#
            prs_temp = imresize3(PRS0,'Scale',[5,5,1],'Method','nearest');
            prs_temp = imresize3(prs_temp,'Scale',[1/11,1/11,1],'Method','box');
            
            if ~isempty(find(prs_temp(:)<0, 1))
                1;
            end
            
            
            PRS2_2 = cat(3,PRS2_2,prs_temp);
            fprintf('Radar_%04d_%02d_%02d finished\n',year,mon,day);
            
        end
        TS = datetime(year,mon,1,0,0,0):1/24:datetime(year,mon,eomday(year,mon),23,0,0);
        
        % SAVE DATA
        % file is currently saved for each year.
        % for 1 year hourly data ([50,50,365*24])
        % <int16> ~ 50 Mb, ok to use after double(*)
        % <double> ~ 175.2 Mb, ok;        
        time = datenum(TS);
        
        % E_temp = imresize(RADAR_info.RegionEE(1,:),'Scale',[1,5],'Method','nearest');
        % N_temp = imresize(RADAR_info.RegionNN(:,1)','Scale',[1,5],'Method','nearest');
        % E = imresize(E_temp,'Scale',[1,1/11],'Method','box')/1000;
        % N = imresize(N_temp,'Scale',[1,1/11],'Method','box')/1000;
        
        
        E = region.minE:2.2:region.minE+2.2*(region.dimE-1);% = RADAR_info.EE(loci,locj);
        N = region.minN:2.2:region.minN+2.2*(region.dimN-1);
        
        pr = permute(PRS2_2,[2,1,3])/scaleF;
        filePath = sfPath;
        fileName = sprintf('pr_nimrod_uk_2.2km_1hr_%04d%02d.nc',year,mon);
        createRadarNetCDF(filePath,fileName,N,E,time,pr)

    end
	
end

end


%% AUXILLARY function

function [PRS,RADAR_info] = getPRS(PRS0,RADAR_info,region,loci,locj)

PRS0 = reshape(PRS0,[RADAR_info.row,RADAR_info.col,size(PRS0,2)]);
PRS = double(PRS0(loci,locj,:));
RADAR_info.RegionEE = RADAR_info.EE(loci,locj);
RADAR_info.RegionNN = RADAR_info.NN(loci,locj);

end

function [loci,locj] = getLoc_IRADAR(RADAR_info,region)

E = [region.minE-1.1, region.minE+2.2*(region.dimE-1)+1.1]*1000;
N = [region.minN-1.1, region.minN+2.2*(region.dimN-1)+1.1]*1000;
totLenE = round(2.2*region.dimE)/1;
totLenN = round(2.2*region.dimN)/1;

rind_E = getRI(E,RADAR_info.E);
rind_N = getRI(N,RADAR_info.N);

% check totLenE
if abs(diff(rind_E))+1 ~= totLenE | abs(diff(rind_N))+1 ~= totLenN
    fprintf('disE in Radar:%02d ',totLenE-abs(diff(rind_E))-1);
    fprintf('disN in Radar:%02d\n',totLenN-abs(diff(rind_N))-1);
    rind_E(2) = rind_E(1)+totLenE-1;
    rind_N(1) = rind_N(2)+totLenN-1;
    % uiwait(msgbox('Chekc totLenE','Title','modal'));
end


loci = rind_N(1):-1:rind_N(2);
locj = rind_E(1):1:rind_E(2);

    function rind = getRI(ran,Vec)
        findi = @(i)find(abs(Vec-ran(i))==min(abs(Vec-ran(i))));
        r1 = findi(1);
        r2 = findi(2);
        rind = [r1(1),r2(end)];
    end
end


