
function [FLO_test,FLO_pred,E,N] = Model_floodPrediction(FOR_test,FOR_analog)

[FLO_test,FLO_pred] = deal(struct('eventNo',[],'floodmaps',[],...
    'unit','1km','mode','maxFilter'));
FLO_test.floodmaps = NaN(0,40,30);%[idx,[Loc]]
FLO_pred.floodmaps = NaN(0,40,30,12);%[idx,[loc],ensemble];


for testInd = 1:length(FOR_test.eventNo)
    
    [FLO_test.floodmaps(testInd,:,:),E,N] = loadFloodMaps(FOR_test.eventNo(testInd),FLO_test.unit,FLO_test.mode);
     
    [FLO_pred.floodmaps(testInd,:,:,:),E,N] = loadFloodMaps(FOR_analog.eventNo(testInd,:),FLO_pred.unit,FLO_pred.mode);   
    
end


FLO_test.eventNo = FOR_test.eventNo;
FLO_pred.eventNo = FOR_analog.eventNo;

end

function [floodmaps,E,N] = loadFloodMaps(eventNo,unit,mode)

arguments
    eventNo (1,:) double
    unit (1,:) char = '1km'
    mode (1,:) char = 'marFilter'
end

floodmaps = NaN(40,30,numel(eventNo));
filePath = 'G:\BIGDATA\TOPIC 2\ProcessedFiles\FloodMaps_40_30';

for evi = 1:numel(eventNo)
    
    fileName = sprintf('FloodMaps_Merged_EventNo%03d.mat',eventNo(evi));
    A = load([filePath,filesep,fileName],'FMap','E','N');
    floodmaps(:,:,evi) = A.FMap;
    
end

E = A.E;
N = A.N;

end

