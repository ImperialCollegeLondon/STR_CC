

function [fea, fea_val, Um, numP] = getImageFeature(DAT,DAT_val,method,pcNum)
% getImageFeature extract feature and reduce the dimnensionality of DAT
%
% Input
%      DAT: 3d # <double> [loc1,loc2,sampleN]
%      DAT: 3d # <double> [loc1,loc2,sampleN]
%      method: <char>
% Output:
%      fea: # <double> [??,??]
%      what else?:
%
% @ Yuting Chen
% Imperial College London
% yuting.chen17@imperial.ac.uk
% Ref: Lu Haiping et al., Uncorrelated Multilinear Principal Component Analysis
%      for Unsupervised Multilinear Subspace Learning
% Required:
% Yuting_Func folder

arguments
    DAT
    DAT_val
    method (1,:) char {mustBeMember(method,{'UMPCA','MPCA','PCA'})}
    pcNum (1,1) double
end

% R1 = load('K:\DATA_FCS\RainEvent\Birm-34-30\EventNo012.mat');
% R2 = load('K:\DATA_FCS\RainEvent\Birm-34-30\EventNo061.mat');
% fea2D = cat(3,R1.RMap,R2.RMap);


switch (method)
    
    case 'UMPCA'
        
        fea2D = DAT;
        N = ndims(fea2D)-1;%Order of the tensor sample
        Is = size(fea2D);%34x30x...
        
        numSpl = Is(3);%There are $numSpl$ image samples
        numP = pcNum;%prod(size(fea2D,[1,2]));
        [Us,TXmean,odrIdx] = UMPCA(fea2D,numP);
        fea2D = fea2D-repmat(TXmean,[ones(1,N), numSpl]);%Centering
        numP = length(odrIdx);
        newfea = zeros(numSpl,numP);
        
        for iP = 1:numP
            projFtr = ttv(tensor(fea2D),Us(:,iP),[1 2]);
            newfea(:,iP) = projFtr.data;
        end
        newfea = newfea(:,odrIdx);% newfea is the final feature vector to be
        %fed into a standard classifier (e.g., nearest neighbor classifier)
        
        fea = newfea;
        Um = Us;
        
        fea2D = DAT_val;
        N = ndims(fea2D)-1;%Order of the tensor sample
        Is = size(fea2D);
        numSpl = Is(3);%There are $numSpl$ image samples
        numP = pcNum;
        fea2D = fea2D-repmat(TXmean,[ones(1,N), numSpl]);%Centering
        numP = length(odrIdx);
        newfea = zeros(numSpl,numP);
        for iP = 1:numP
            projFtr = ttv(tensor(fea2D),Us(:,iP),[1 2]);
            fea_val(:,iP) = projFtr.data;
        end
        fea_val = fea_val(:,odrIdx);
        fea_val = fea_val(:,1:pcNum);
        
    case 'PCA'
        
        preTrans = @(x)reshape(log(x+0.01),numel(x(:,:,1)),[]);% Box-Cox Transformation
        XTrain = transpose(preTrans(DAT));
        % Preprocessing of PCA func:
        % original ingredients data centered by subtracting the column means from corresponding columns.
        [coeff,scoreTrain,latent,tsquared,explained,mu]  = pca(XTrain,'NumComponents',pcNum);
        U_pm = coeff(:,1:pcNum);% at least 95% variance is preserved
        
        X_lp = transpose(preTrans(DAT))-mu;% exp(transpose(preTrans(RMapS)))-0.01;
        fea = X_lp*U_pm;
        Um = U_pm;
        numP = pcNum;
        
        X_lp = transpose(preTrans(DAT_val))-mu;% exp(transpose(preTrans(RMap)))-0.01;
        fea_val = X_lp*U_pm;
        
        
    otherwise
        
end
end
