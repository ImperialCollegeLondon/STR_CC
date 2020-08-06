ENSEMBLENO = getEnsNos();

REGIONS = REGIONS_info();
IntFac = 32; % to make the resolution = (Nimrod)
options = 'Save';
for ENSEMBLENO_1 = ENSEMBLENO
    for MON = 6:8
        
        for region = [REGIONS.CPM_NE,REGIONS.CPM_NW,REGIONS.CPM_S]
            % [REGIONS.SCO,REGIONS.WAL,REGIONS.EUK]
            % [REGIONS.SWestuk,REGIONS.Westuk,REGIONS.Scotland,REGIONS.London]%
            
            data = struct('Years',[1980,2000],...
                'fileGetPath','K:/UkCp18/',...
                'savePath',['D:/UKCP18/',region.Name]);
            readCPM_pr(region,ENSEMBLENO_1,MON,IntFac,data,options);
            
            %             data = struct('Years',[2020,2040],...
            %                 'fileGetPath','K:/UkCp18_FutureTemp/',...
            %                 'savePath',['D:/UKCP18_Future/',region.Name]);
            %             readCPM_pr(region,ENSEMBLENO_1,MON,IntFac,data,options);
            
            data = struct('Years',[2060,2080],...
                'fileGetPath','K:/UkCp18_FutureTemp/',...
                'savePath',['D:/UKCP18_Future_2060_2080/',region.Name]);
            readCPM_pr(region,ENSEMBLENO_1,MON,IntFac,data,options);
            
        end
    end
end

%%
fprintf('EXTRACTION Done!\n')
ProcessIDAFs_AMA
% UKCP_process_ConvStorm_JJA