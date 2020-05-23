% This script is to see the effect of domain size on extracted IMF series
% and WAR series.

clear;clc
[T1,IMF1,WAR1]= getIt(6);
[T2,IMF2,WAR2]= getIt(7);
[T3,IMF3,WAR3]= getIt(8);

T = [T1;T2;T3];
IMF = [IMF1;IMF2;IMF3];
WAR = [WAR1;WAR2;WAR3];
clear T1 T2 T3 IMF1 IMF2 IMF3 WAR1 WAR2 WAR3

% save('STATS_WALES_2060-2080_Ens05.mat','T','IMF','WAR');
figure;
subplot 211
plot(T(:,1),IMF(:,1));ylim([0,2])
subplot 212
hold on
plot(T(:,1),WAR(:,1));ylim([0,1])
ax = gca;
plot(ax.XLim,[0.05,0.05]);
title('BIGER AREA(*00Km**00Km)')

[T1,IMF1,WAR1]= getsmall(6);
[T2,IMF2,WAR2]= getsmall(7);
[T3,IMF3,WAR3]= getsmall(8);

T = [T1;T2;T3];
IMF = [IMF1;IMF2;IMF3];
WAR = [WAR1;WAR2;WAR3];
clear T1 T2 T3 IMF1 IMF2 IMF3 WAR1 WAR2 WAR3

% save('STATS_SMALLWALES_2060-2080_Ens05.mat','T','IMF','WAR');

figure;
subplot 211
plot(T(:,1),IMF(:,1));ylim([0,2])
subplot 212
hold on
plot(T(:,1),WAR(:,1));ylim([0,1])
ax = gca;
plot(ax.XLim,[0.02,0.02]);
title('SMALLER AREA(110Km*110Km)')


function [T,IMF,WAR]= getsmall(mon)
load(['D:\UKCP18_Future_2060_2080\Westuk\',sprintf('Ensems_mon%02d.mat',mon)])
A = (RainEnsembles{3});
IMF = (squeeze(nansum(nansum(double(A),1),2)))/32./prod(size(A,[1,2]));
WAR = (squeeze(nansum(nansum(double(A)>0.1*32,1),2)))./prod(size(A,[1,2]));
[hh,dd,yy] = meshgrid(1:24,1:30,2061:2080);
T = datetime(yy,mon,dd,hh,0,0);
T = permute(T,[2,1,3]);
T = T(:);
reshapeIt = @(x)reshape(x,24*30,20);
IMF = reshapeIt(IMF);
T = reshapeIt(T);
WAR = reshapeIt(WAR);
end

function [T,IMF,WAR]= getIt(mon)
load(['D:\UKCP18_Future_2060_2080\bigWAL\',sprintf('Ensems_05_mon%02d.mat',mon)])
A = (RainEnsembles{1, 1});
IMF = (squeeze(nansum(nansum(double(A),1),2)))/32./prod(size(A,[1,2]));
WAR = (squeeze(nansum(nansum(double(A)>0.1*32,1),2)))./prod(size(A,[1,2]));
[hh,dd,yy] = meshgrid(1:24,1:30,2061:2080);
T = datetime(yy,mon,dd,hh,0,0);
T = permute(T,[2,1,3]);
T = T(:);
reshapeIt = @(x)reshape(x,24*30,20);
IMF = reshapeIt(IMF);
T = reshapeIt(T);
WAR = reshapeIt(WAR);
end