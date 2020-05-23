function [version,testConfig,filefolder] = getExpNoInfo(expNo)
% getExpNoInfo(expNo) will get Corresponding exp info,
% see detail in report WP3.3

arguments
    expNo (1,:) char
end
testConfig = struct();
testConfig.IDW = false;% false is better (at least for event 90)
testConfig.uniqueEvent = false;
testConfig.includePrevious = false;

switch expNo
    case 'C1P5-1'
        version = 'v1';
        testConfig.imageMethod = 'PCA';%'UMPCA';% {mustBeMember(method,{'UMPCA','MPCA','PCA'})}
        testConfig.pcNum = 5;
        testConfig.rainfallSource = 'KED';
        filefolder = ['K:\DATA_FCS\CrossVal\',sprintf('%s_%s_PCs%02d%s',version,testConfig.imageMethod,testConfig.pcNum,filesep)];
    case 'C1P5-2'
        version = 'v1';
        testConfig.imageMethod = 'PCA';%'UMPCA';% {mustBeMember(method,{'UMPCA','MPCA','PCA'})}
        testConfig.pcNum = 5;
        testConfig.rainfallSource = 'KED';
        testConfig.excludeWAR = true;
        filefolder = ['K:\DATA_FCS\CrossVal\',sprintf('%s_%s_PCs%02d_noWAR%s',version,testConfig.imageMethod,testConfig.pcNum,filesep)];
    case 'C1U5-1'
        version = 'v1';
        testConfig.imageMethod = 'UMPCA';%'UMPCA';% {mustBeMember(method,{'UMPCA','MPCA','PCA'})}
        testConfig.pcNum = 5;
        testConfig.rainfallSource = 'KED';
        filefolder = ['K:\DATA_FCS\CrossVal\',sprintf('%s_%s_PCs%02d%s',version,testConfig.imageMethod,testConfig.pcNum,filesep)];
        
    case 'C1U10-1'
        version = 'v1';
        testConfig.imageMethod = 'UMPCA';%'UMPCA';% {mustBeMember(method,{'UMPCA','MPCA','PCA'})}
        testConfig.pcNum = 10;
        testConfig.rainfallSource = 'KED';
        filefolder = ['K:\DATA_FCS\CrossVal\',sprintf('%s_%s_PCs%02d%s',version,testConfig.imageMethod,testConfig.pcNum,filesep)];
        
    case 'C2P5-Rad1'
        version = 'v2';
        testConfig.imageMethod = 'PCA';%'UMPCA';% {mustBeMember(method,{'UMPCA','MPCA','PCA'})}
        testConfig.pcNum = 5;
        testConfig.rainfallSource = 'RAD';
        filefolder = ['K:\DATA_FCS\CrossVal\',sprintf('%s_%s_PCs%02d_Rain%s%s',...
            version,testConfig.imageMethod,testConfig.pcNum,...
            testConfig.rainfallSource,filesep)];
        
    case 'C1P5In-1'
        version = 'v1';
        testConfig.imageMethod = 'PCA';%'UMPCA';% {mustBeMember(method,{'UMPCA','MPCA','PCA'})}
        testConfig.pcNum = 5;
        testConfig.rainfallSource = 'KED';
        testConfig.analogtimegap = 3;
        filefolder = ['K:\DATA_FCS\CrossVal\',sprintf('%s_%s_PCs%02d_tGap%01d%s',...
            version,testConfig.imageMethod,testConfig.pcNum,testConfig.analogtimegap,filesep)];
        
    case 'C1P5-Rad1'
        version = 'v1';
        testConfig.imageMethod = 'PCA';%'UMPCA';% {mustBeMember(method,{'UMPCA','MPCA','PCA'})}
        testConfig.pcNum = 5;
        testConfig.rainfallSource = 'RAD';
        filefolder = ['K:\DATA_FCS\CrossVal\',sprintf('%s_%s_PCs%02d_Rain%s%s',...
            version,testConfig.imageMethod,testConfig.pcNum,...
            testConfig.rainfallSource,filesep)];
        
    case 'C1P5In-Rad1'
        version = 'v1';
        testConfig.imageMethod = 'PCA';%'UMPCA';% {mustBeMember(method,{'UMPCA','MPCA','PCA'})}
        testConfig.pcNum = 5;
        testConfig.rainfallSource = 'RAD';
        testConfig.analogtimegap = 3;
        filefolder = ['K:\DATA_FCS\CrossVal\',sprintf('%s_%s_PCs%02d_Rain%s_tGap%01d%s',...
            version,testConfig.imageMethod,testConfig.pcNum,...
            testConfig.rainfallSource,testConfig.analogtimegap,filesep)];
    case 'C1U5-Rad1'
        version = 'v1';
        testConfig.imageMethod = 'UMPCA';%'UMPCA';% {mustBeMember(method,{'UMPCA','MPCA','PCA'})}
        testConfig.pcNum = 5;
        testConfig.rainfallSource = 'RAD';
        filefolder = ['K:\DATA_FCS\CrossVal\',sprintf('%s_%s_PCs%02d_Rain%s%s',...
            version,testConfig.imageMethod,testConfig.pcNum,...
            testConfig.rainfallSource,filesep)];
    case 'C1U10-Rad1'
        version = 'v1';
        testConfig.imageMethod = 'UMPCA';%'UMPCA';% {mustBeMember(method,{'UMPCA','MPCA','PCA'})}
        testConfig.pcNum = 10;
        testConfig.rainfallSource = 'RAD';
        filefolder = ['K:\DATA_FCS\CrossVal\',sprintf('%s_%s_PCs%02d_Rain%s%s',...
            version,testConfig.imageMethod,testConfig.pcNum,...
            testConfig.rainfallSource,filesep)];
    otherwise
        error('This ExpNo doesnt exist!\n');
end

end