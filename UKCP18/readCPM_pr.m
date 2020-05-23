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
    MON (1,1) double
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
[region.i,region.j] = getRegionIJ(E,N,region.minE,region.minN);

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
        
        if strcmpi(options,'Mean')
            MRain = [];
        else
            Rain = NaN(region.dimE,region.dimN,20*30*24);
        end
        
        
        L = 1;
        
        for year = data.Years(1):data.Years(end)
            tic
            if ~((year==data.Years(1) && mon ~=12) || (year == data.Years(end) && mon == 12))
                
                fileName = sprintf('pr_rcp85_land-cpm_uk_2.2km_%s_1hr_%04d%02d01-%04d%02d30.nc',...
                    ENSEMBLENO{M},year,mon,year,mon);
                
                listaRain = dir([filePath,fileName]);
                
                
                if ~isempty(listaRain)
                    
                    listaRain = fullfile({listaRain.folder},{listaRain.name});
                    listaRain = listaRain{1};
                    A = ncinfo(listaRain);
                    
                    if strcmpi(options,'Mean')
                        MRain(:,:,L) = nanmean(squeeze(ncread(listaRain,'pr',...
                            [region.i,region.j,1,1],...
                            [region.dimE,region.dimN,...
                            A.Variables(1).Dimensions(:,3).Length,...
                            A.Variables(1).Dimensions(:,4).Length])),3);
                    else
                        Rain(:,:,(L-1)*30*24+1:(L)*30*24) = squeeze(ncread(listaRain,'pr',...
                            [region.i,region.j,1,1],...
                            [region.dimE,region.dimN,...
                            A.Variables(1).Dimensions(:,3).Length,...
                            A.Variables(1).Dimensions(:,4).Length]));
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
            if getfield(whos('Rain'),'bytes') > 5e9
                uiwait(msgbox('Almost out of memory, hope next time you can seperate data into files',...
                    'Oooops','__'));
            end
            RainEnsembles{M} = Rain;
            fprintf('Ensembles %02d\n',M);
        end
        
    end
    
    if strcmpi(options,'save')
        RainEnsembles = cellfun(@(x) int16(round(IntFac*x)), RainEnsembles, ...
            'UniformOutput',false);
        save([saveDir,'.mat'],'RainEnsembles','-v7.3');
    else
        IntFac = 1;
        fprintf('Mon %02d Exported out, in <double> format\n',M);
    end
    
    fprintf('---%s-Mon%02d---\n',region.Name,mon);
end

end