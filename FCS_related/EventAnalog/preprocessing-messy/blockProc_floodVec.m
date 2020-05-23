
FP = dir('G:\BIGDATA\TOPIC 2\ProcessedFiles\WCMaps_accurate\WCMaps*.mat');
SP = 'G:\BIGDATA\TOPIC 2\ProcessedFiles\WCMaps\';

mkdir('G:\BIGDATA\TOPIC 2\ProcessedFiles\WCMaps_100Max\')

for fil = FP'
    
    filename = [fil.folder,filesep,fil.name];
    WCM = load(filename,'FloodVec');
    
    
    FloodVec = WCM.FloodVec(1:64);
    save(['G:\BIGDATA\TOPIC 2\ProcessedFiles\WCMaps_100Max\',fil.name],'FloodVec');
    
%     blocksize = [100,1];
%     FloodVec = blockproc(WCM.FloodVec, blocksize, @(x)max(x.data(:)));
%     
%     save(['G:\BIGDATA\TOPIC 2\ProcessedFiles\WCMaps_100Max\',fil.name],'FloodVec');
    
end



% load('G:\BIGDATA\TOPIC 2\ProcessedFiles\WCMaps\WCMaps_Event051_WCNo11.mat')

%%


FP = dir('G:\BIGDATA\TOPIC 2\ProcessedFiles\WCMaps_accurate\WCMaps*_WCNo11.mat');
SP = 'G:\BIGDATA\TOPIC 2\ProcessedFiles\WCMaps\';

mkdir('G:\BIGDATA\TOPIC 2\ProcessedFiles\WCMaps_100Max\')

for fil = FP'
    
    filename = [fil.folder,filesep,fil.name];
    WCM = load(filename,'FloodVec');
    
    
    FloodVec = WCM.FloodVec(1:64);
    plot(FloodVec);
    hold on;
    % save(['G:\BIGDATA\TOPIC 2\ProcessedFiles\WCMaps_100Max\',fil.name],'FloodVec');
    
%     blocksize = [100,1];
%     FloodVec = blockproc(WCM.FloodVec, blocksize, @(x)max(x.data(:)));
%     
%     save(['G:\BIGDATA\TOPIC 2\ProcessedFiles\WCMaps_100Max\',fil.name],'FloodVec');
    
end


