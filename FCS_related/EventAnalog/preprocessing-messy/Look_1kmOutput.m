
FP = dir('G:\BIGDATA\TOPIC 2\ProcessedFiles\FloodMaps\FloodMaps*.mat');
SP = 'G:\BIGDATA\TOPIC 2\ProcessedFiles\WCMaps\';

mkdir('G:\BIGDATA\TOPIC 2\ProcessedFiles\WCMaps_1km_Max\')


R = zeros(40,30);
for fil = FP'
    
    filename = [fil.folder,filesep,fil.name];
    WCM = load(filename,'FMap','E','N');
    
    
    
    FMap = WCM.FMap;
    % save(['G:\BIGDATA\TOPIC 2\ProcessedFiles\WCMaps_1km_Max\',fil.name],'FMap');
    R = R+FMap;
    
    
    pcolor(WCM.E,WCM.N,WCM.FMap);shading flat;
    cptcmap('GMT_no_green', 'mapping','scaled');caxis([0,100])
    title(upper(fil.name(end-13:end-4)))
    colorbar
    xlabel('Easting')
    ylabel('Northing')
    % pause(0.01)
    
    saveas(gcf,sprintf('C:/Users/Yuting Chen/Dropbox (Personal)/Data_PP/Fig_FCS/output_1km_FMaps/%s',fil.name(1:end-4)),'png')


end
%%
R(R==0)=NaN;
pcolor(R)
shading flat
cptcmap('diff_darkBlue_darkRed', 'mapping','scaled');