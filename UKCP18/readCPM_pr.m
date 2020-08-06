%--------------------------------------------------%
% this file is to read ensembles of CPM Pr.
%
% @ Yuting Chen
% Imperial College London
% yuting.chen17@imperial.ac.uk
% Update: 2019.12.**
%--------------------------------------------------%


function [RainEnsembles,IntFac] = readCPM_pr(region,ENSEMBLENO,MON,IntFac,data,options)
% READCPM_PR() gives .....
%
% Input:region
%       ENSEMBLENO: (1,>=1) <cell>
%       MON:
%       IntFac:recommend: 32 (format of NIMROD RADAR)
%       data:
%       options:
% Output:RainEnsembles
%        IntFac
%
% Example 1: save to file
%
%         ENSEMBLENO=getEnsNos();
%         IntFac = 32; % to make the resolution = (Nimrod)
%         REGIONS = REGIONS_info();
%         options = 'Save';
%
%         region = REGIONS.SWestuk;
%         data = struct('Years',[2020,2040],...
%             'fileGetPath','K:/UkCp18_FutureTemp/',...
%             'savePath',['D:/UKCP18_Future/',region.Name]);
%         readCPM_pr(REGIONS.SWestuk,ENSEMBLENO,MON,IntFac,data,options);
%
% Example 2: output as mean
%         clear;clc
%         IntFac = 1;
%         options = 'Mean';
%         [RainEnsembles,IntFac] = readCPM_pr(region,ENSEMBLENO,MON,IntFac,data,options);
%
%
% # Currently available forma for $data$ <struct> #
% data = struct('Years',[1980,2000],...
%             'fileGetPath','K:/UkCp18/',...
%             'savePath',['D:/UKCP18/',region.Name]);
% data = struct('Years',[2020,2040],...
%             'fileGetPath','K:/UkCp18_FutureTemp/',...
%             'savePath',['D:/UKCP18_Future/',region.Name]);
% data = struct('Years',[2060,2080],...
%             'fileGetPath','K:/UkCp18_FutureTemp/',...
%             'savePath',['D:/UKCP18_Future_2060_2080/',region.Name]);
%
%
% @ Yuting Chen
% Imperial College London


arguments
    
    region (1,1) struct
    ENSEMBLENO (1,:) cell
    MON (1,:) double
    IntFac (1,1) double = 32;
    data (1,1) struct = struct('Years',[1980,2000],...
        'fileGetPath','K:/UkCp18/',...
        'savePath',['D:/UKCP18/',region.Name]);
    options (1,:) char {mustBeMember(options,{'Save','Mean'})} = 'save'
    
end


%% Configuration
filePath = ['K:/UkCp18/',ENSEMBLENO{1},'/pr_rcp85_land-cpm_uk_2.2km_'];
ncFileName = [filePath,ENSEMBLENO{1},'_1hr_19910801-19910830.nc'];
LAT=ncread(ncFileName,'latitude');
LON=ncread(ncFileName,'longitude');
[E, N] = ll2os(LAT, LON);
E = E/1000;
N = N/1000;

mkdir(data.savePath);
fprintf(sprintf('%04d-%04d data will be extracted',data.Years(1),data.Years(end)));

%%
if isRegularRegion(region)
    [region.i,region.j] = getRegionIJ(E,N,region.minE,region.minN);
    if region.dx ~= 2.2
        region.dimE = round(region.dimE*region.dx/2.2);
        region.dimN = round(region.dimN*region.dx/2.2);
    end
    i_this = region.i:region.i+region.dimE-1;
    j_this = region.j:region.j+region.dimN-1;
    Rain_E = E(i_this,j_this);
    Rain_N = N(i_this,j_this);
else
    % region is represented by a polygon {region.E},{region.N}
    % (rectangular shape having an rotation angle ~= 0)
    angle = atan((region.E(2)-region.E(1))/(region.N(2)-region.N(1)))*180/pi;
    E = imrotate(E,angle,'nearest','loose');
    N = imrotate(N,angle,'nearest','loose');
    [region.i,region.j] = arrayfun(@(e0,n0)getRegionIJ(E,N,e0,n0),region.E,region.N);
    i_this = min(region.i):max(region.i);
    j_this = min(region.j):max(region.j);
    Rain_E = E(i_this,j_this);
    Rain_N = N(i_this,j_this);
end

warning off

for mon = MON
    
    RainEnsembles = cell(length(ENSEMBLENO),1);
    
    if length(ENSEMBLENO)>1
        saveDir = [data.savePath,'/Ensems_mon',sprintf('%02d',mon)];
    else
        saveDir = [data.savePath,'/',sprintf('Ensems_%s_mon%02d',ENSEMBLENO{1},MON)];
    end
    for M=1:length(ENSEMBLENO)
        
        filePath = [data.fileGetPath,ENSEMBLENO{M},'/'];
        if isRegularRegion(region)
            if strcmpi(options,'Mean')
                MRain = [];
            else
                Rain = NaN(region.dimE,region.dimN,20*30*24);
            end
        else
            if strcmpi(options,'Mean')
                MRain = [];
            else
                Rain = [];
            end
        end
        
        L = 1;
        
        for year = data.Years(1):data.Years(end)
            tic
            if ~((year==data.Years(1) && mon ~=12) || (year == data.Years(end) && mon == 12))
                
                fileName = sprintf('pr_rcp85_land-cpm_uk_2.2km_%s_1hr_%04d%02d01-%04d%02d30.nc',...
                    ENSEMBLENO{M},year,mon,year,mon);
                
                thisNCFile = dir([filePath,fileName]);
                
                if ~isempty(thisNCFile)
                    
                    thisNCFile = fullfile({thisNCFile.folder},{thisNCFile.name});
                    thisNCFile = thisNCFile{1};
                    A = ncinfo(thisNCFile);
                    
                    if isRegularRegion(region)
                        Rain_this = squeeze(ncread(thisNCFile,'pr',...
                            [region.i,region.j,1,1],...
                            [region.dimE,region.dimN,...
                            A.Variables(1).Dimensions(:,3).Length,...
                            A.Variables(1).Dimensions(:,4).Length]));
                        if strcmpi(options,'Mean')
                            MRain(:,:,L) = nanmean(Rain_this,3);
                        else
                            Rain(:,:,(L-1)*30*24+1:(L)*30*24) = Rain_this;
                        end
                    else
                        Rain_this = squeeze(ncread(thisNCFile,'pr'));
                        Rain_this = imrotate(Rain_this,angle,'nearest','loose');
                        Rain_this = Rain_this(i_this,j_this,:);
                        if strcmpi(options,'Mean')
                            MRain(:,:,L) = nanmean(Rain_this,3);
                        else
                            Rain = cat(3,Rain,Rain_this);
                        end
                    end
                else
                    
                end
                L = L+1;
                
            end
            toc
        end
        
        if strcmpi(options,'UKMean')
            RainEnsembles{M} = MRain;
            fprintf('Ensembles %02d\n',M);
        else
            if getfield(whos('Rain'),'bytes') > 8e9
                uiwait(msgbox('@YC: Almost out of memory! Hope next time you can seperate data and save into several files',...
                    'Oooops','__'));
            end
            RainEnsembles{M} = Rain;
            fprintf('Ensembles %02d\n',M);
        end
        
    end
    
    if strcmpi(options,'save')
        RainEnsembles = cellfun(@(x) int16(round(IntFac*x)), RainEnsembles, ...
            'UniformOutput',false);
        save([saveDir,'.mat'],'RainEnsembles','Rain_E','Rain_N','-v7.3');
    else
        IntFac = 1;
        fprintf('Mon %02d Exported out, in <double> format\n',M);
    end
    
    fprintf('---%s-Mon%02d---\n',region.Name,mon);
end

end