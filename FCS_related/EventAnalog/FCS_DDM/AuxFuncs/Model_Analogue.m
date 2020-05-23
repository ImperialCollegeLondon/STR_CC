function [RES_analog,DIS_analogL1,DIS_analogL2] = Model_Analogue(AnalogInput,FANum)
% Model_Analogue will do two-layer analog given AnalogInput and FANum
%
% Input: AnalogInput:
%        FANum:
% Output:RES_analog:
%        DIS_analogL1:
%        DIS_analogL2:
%
% @ Yuting Chen
% yuting.chen17@imperial.ac.uk

% getinput
PRED = AnalogInput.PRED;
PREDval = AnalogInput.PREDval;
PREDRAIN = AnalogInput.PREDRAIN;
PREDRAINval = AnalogInput.PREDRAINval;
raintime = AnalogInput.raintime;
raintimeval = AnalogInput.raintimeval;
timegap = AnalogInput.timegap; % unit [1h]

% FIRST LAYER ----- FORCING ANALOGUES
FA = PRED;
FAval = PREDval;


[coeff,scoreTrain,latent,tsquared,explained,mu] = pca(FA,'NumComponents',size(FA,2));
U_pm = coeff(:,1:size(FA,2));% 95% variance is preserved
FA = FA*U_pm;
FAval = FAval*U_pm;


[Idx_1, Dist_1_all] = firstLayer(FA,FAval,FANum,raintime,timegap);

% SECOND LAYER ----- RAINFALL ANALOGUES (Based on PCA)

Idx_2 = [];

FA = PREDRAIN;
FAval = PREDRAINval;

Dist_2 = [];

for testi = 1:size(Idx_1,1)
    
    FA_temp = FA(Idx_1(testi,:),:);
    FAval_temp = FAval(testi,:);
    
    [Idx_2(testi,:),Dist_2(testi,:)] = secondLayer(FA_temp,FAval_temp,12);
    
    Dist_1(testi,:) = Dist_1_all(testi,Idx_2(testi,:));
    [Idx_2(testi,:)] = Idx_1(testi,Idx_2(testi,:));% store the original id.
    
end

% only save necessary output
RES_analog = Idx_2;
DIS_analogL1 = Dist_1;
DIS_analogL2 = Dist_2;

end

function  [IdxES, IdxDist] = firstLayer(FA,FAval,FANum,raintime,timegap)

rng('shuffle')
% n = size(FA,1);
% idx = randsample(n,500);
% X = FA(~ismember(1:n,idx),:); % Training data
% Y = FA(idx,:); % Test data

X = FA; % Training dataset
Y = FAval; % Test data
MdlES = ExhaustiveSearcher(X);

if timegap == 0
    
    [IdxES,IdxDist] = knnsearch(MdlES,Y,'K',FANum,'Distance','seuclidean');

else
    
    [IdxES,IdxDist] = knnsearch(MdlES,Y,'K',min(FANum*timegap*(60/5)*2,size(X,1)),...
        'Distance','seuclidean');
    
    % Dependent candidates will be all excluded.
    for testi = 1:size(Y,1)
        for faNo =  1:FANum % repeat for FANum times
            
            excInd = find(raintime(IdxES(testi,:)) > (raintime(IdxES(testi,faNo))-timegap/24) & ...
                raintime(IdxES(testi,:)) < (raintime(IdxES(testi,faNo))+timegap/24) & ...
                raintime(IdxES(testi,:)) ~= (raintime(IdxES(testi,faNo))));
            
            if ~isempty(excInd) && ~isempty(excInd(excInd>faNo))
                % set it to be extremely large instead of changing the the
                % vector size.
                % in order to exclude it in the final step.
                excInd = excInd(excInd>faNo);
                IdxDist(testi,excInd) = IdxDist(testi,excInd)*0+1e6;
                
            end
        end
    end
    [IdxDist,I] = sort(IdxDist,2);
    IdxES = IdxES(I);
    IdxDist = IdxDist(:,1:FANum);
    IdxES = IdxES(:,1:FANum);
end
end

function  [IdxES,IdxDist] = secondLayer(FA,FAval,FANum)

rng('shuffle')
% n = size(FA,1);
% idx = randsample(n,500);
% X = FA(~ismember(1:n,idx),:); % Training data
% Y = FA(idx,:); % Test data

X = FA; % Training dataset
Y = FAval; % Test data

% MdlKDT = KDTreeSearcher(X);
% IdxKDT = knnsearch(MdlKDT,Y,'K',FANum);%faster
MdlES = ExhaustiveSearcher(X);

[IdxES,IdxDist] = knnsearch(MdlES,Y,'K',FANum,'Distance','seuclidean');

end
