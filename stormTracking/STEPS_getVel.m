% DATA FORMAT of vel_x,vel_y,R:
% 2d # <double>
% [northing, easting]
% left-lower CORNER of each # represents to:
%            [min northing, min easting]


vel_x = load('H:\CODE_PYTHON\temp_try\test_vel_x.csv');
vel_y = load('H:\CODE_PYTHON\temp_try\test_vel_y.csv');
R = load('H:\CODE_PYTHON\temp_try\test_R0.csv');

load('D:\UKCP18\bigSCO\oneM.mat','R')

R(R==0) = NaN;
for tim = 1:1:19
subplot(1,2,1)
pcolor(squeeze(R(tim,:,:)));shading flat
cptcmap('precip_meteoswiss', 'mapping','direct');%,'ncol',20);
set(gca,'YDir','Reverse')
xlabel('j index');ylabel('i index')

subplot(1,2,2)
quiver(imresize(vel_x,1/25,'bilinear'),...
    imresize(vel_y,1/25,'bilinear'))

set(gca,'YDir','Reverse')
xlabel('j index');ylabel('i index')

pause(0.1);

end