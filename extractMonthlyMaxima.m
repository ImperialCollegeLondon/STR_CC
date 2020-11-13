clear;clc
close all
%%
% data = struct('Years',[1980,2000],...
%     'fileGetPath','K:/UkCp18/',...
%     'savePath',['D:/UKCP18/','am']);
% options = 'Max';
% readCPM_pr_annualMax(getEnsNos(),6:8,32,data,options);
%%
load('D:\UKCP18\am\Ensems.mat');
[parmhat] = deal(cell(1,numel(RainEnsembles)));
for ensNo = 1:numel(RainEnsembles)
    parmhat{ensNo} = NaN(prod(size(RainEnsembles{ensNo},[1,2])),3);
end
parfor ensNo = 1:numel(RainEnsembles)
    tic
    Rain = RainEnsembles{ensNo};Rain = reshape(Rain,[prod(size(Rain,[1,2])),size(Rain,3)]);
    for i = 1:size(Rain,1)
        ts = Rain(i,:);
        [parmhat{ensNo}(i,:),hobs,hfit] = fitGEV(ts,'method','Gringorten');
        % I{ensNo}(i,1) = gevinv(P,parmhat{ensNo}(i,1),parmhat{ensNo}(i,2),parmhat{ensNo}(i,3));
    end
    % I{ensNo} = reshape(I{ensNo},size(RainEnsembles{1},[1,2]));
    toc
    ensNo
end
save('D:\UKCP18\UK\AM_param_CPM_1980-2000.mat','parmhat');
%%
load('D:\UKCP18\am\Ensems.mat');
load('D:\UKCP18\UK\AM_param_CPM_1980-2000.mat','parmhat');
for RT = [2,5,10]
    tic
    P = 1-(1./RT);
    [I] = deal(cell(1,numel(RainEnsembles)));
    for ensNo = 1:numel(RainEnsembles)
        I{ensNo} = NaN(prod(size(RainEnsembles{ensNo},[1,2])),1);
    end
    parfor ensNo = 1:numel(RainEnsembles)
        for i = 1:size(I{ensNo},1)
            I{ensNo}(i,1) = gevinv(P,parmhat{ensNo}(i,1),parmhat{ensNo}(i,2),parmhat{ensNo}(i,3));
        end
        I{ensNo} = reshape(I{ensNo},size(RainEnsembles{1},[1,2]));
    end
    fileName = ['D:\UKCP18\UK\AM_CPM_1980-2000_RT',num2str(RT),'.mat'];
    var = 'I';
    str = 'save(fileName,var);';
    eval(str);
    toc
end
%%
