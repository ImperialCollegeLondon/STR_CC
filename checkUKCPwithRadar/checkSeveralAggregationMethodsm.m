% CHECK result from IMRESIZE3 and BUFFER MEAN

clear;clc
x1 = 2175;x2 = 1725;
data = int16(exp(3*rand(x1,x2,288))*32);
data0 = data;

%% reshape method

data = reshape(data,x1*x2,288);

tic
data = reshape(data,x1*x2,288);
AD2 = 1/12*double(squeeze(nansum(reshape(data,size(data,1),12,[]),2)));

toc
AD2 = reshape(AD2,x1,x2,[]);



%% buffer mean
% This is time consuming for BIG Matrix.
% data = reshape(data,x1*x2,288);
% 
% tic
% AD = [];
% for ii = 1:size(data,1)
%     AD(ii,:) = reshape(nanmean(buffer(double(data(ii,:)),12),1),1,[]);
% end
% toc
% AD = reshape(AD,x1,x2,[]);

%%
tic
AD0 = imresize(double(data),'Scale',[1,1/12]);
toc
AD0 = reshape(AD0,x1,x2,[]);

%%
difA = AD-AD0;
xx = squeeze(AD(1,1,:));
yy = squeeze(AD0(1,1,:));
zz = squeeze(AD2(1,1,:));
plot(xx);hold on; plot(yy);hold on;
plot(zz);


