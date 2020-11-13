%--------------------------------------------------%
% this file is to read ensembles of CPM Pr annual maxima
%
% @ Titian Chen
% Imperial College London
% Update: 2019.12.**
%--------------------------------------------------%


function [RainEnsembles,IntFac] = readCPM_pr_annualMax(ENSEMBLENO,MON,IntFac,data,options)
% READCPM_PR() gives .....
%
% Input:ENSEMBLENO: (1,>=1) <cell>
%       MON:
%       IntFac:recommend: 32 (format of NIMROD RADAR)
%       data:
%       options:
% Output:RainEnsembles
%        IntFac
%
% Example:
% data = struct('Years',[1980,2000],...
%     'fileGetPath','K:/UkCp18/',...
%     'savePath',['D:/UKCP18/','am']);
% options = 'Max';
% readCPM_pr_annualMax(getEnsNos(),6:8,32,data,options);
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
    
    ENSEMBLENO (1,:) cell
    MON (1,:) double
    IntFac (1,1) double = 32;
    data (1,1) struct = struct('Years',[1980,2000],...
        'fileGetPath','K:/UkCp18/');
    options (1,:) char {mustBeMember(options,{'Save','Max'})} = 'save'
    
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

warning off

RainEnsembles = cell(length(ENSEMBLENO),1);

if length(ENSEMBLENO)>1
    saveDir = [data.savePath,'/Ensems'];
else
    saveDir = [data.savePath,'/',sprintf('Ensems_%s',ENSEMBLENO{1})];
end
for M=1:length(ENSEMBLENO)
    
    filePath = [data.fileGetPath,ENSEMBLENO{M},'/'];
    if strcmpi(options,'Max')
        MRain = [];
    end
    
    L = 1;
    
    for year = data.Years(1):data.Years(end)
        itag = 1; thisMonthMax = [];
        tic
        for mon = MON
            if ~((year==data.Years(1) && mon ~=12) || (year == data.Years(end) && mon == 12))
                
                fileName = sprintf('pr_rcp85_land-cpm_uk_2.2km_%s_1hr_%04d%02d01-%04d%02d30.nc',...
                    ENSEMBLENO{M},year,mon,year,mon);
                
                thisNCFile = dir([filePath,fileName]);
                
                if ~isempty(thisNCFile)
                    thisNCFile = fullfile({thisNCFile.folder},{thisNCFile.name});
                    thisNCFile = thisNCFile{1};
                    A = ncinfo(thisNCFile);
                    Rain_this = squeeze(ncread(thisNCFile,'pr',...
                        [1,1,1,1],[inf,inf,inf,inf]));
                    if strcmpi(options,'Max')
                        thisMonthMax(:,:,itag) = nanmax(Rain_this,[],3);
                        itag = itag+1;
                    end     
                end 
            end
        end
        if itag>1 && strcmpi(options,'Max')
            MRain(:,:,L) = nanmax(thisMonthMax,[],3);
            L = L+1;
        end
        tic
    end
    
    if strcmpi(options,'Max')
        RainEnsembles{M} = MRain;
        fprintf('Ensembles %02d\n',M);
        save([saveDir,'.mat'],'RainEnsembles');
    end
end

end