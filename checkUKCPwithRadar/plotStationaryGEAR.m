
source = ['K:\GEAR\CEH_GEAR_daily_GB_',num2str(2010),'.nc'];
finfo = ncinfo(source);
varname = 'x';
x = ncread(source,varname); % easting-OSGB36 Grid reference
varname = 'y';
y = ncread(source,varname); % northing-OSGB36 Grid reference

load('StationaryCheck_UK_bat_1_386.mat')
RES = [];
tag = 1;
for bat = 1:386
    try
    Y_coor = y(3*bat+1 : 3*bat+3);
    X_coor = x(1:end);
    [X_coor,Y_coor] = meshgrid(X_coor,Y_coor);
    XX = X_coor;
    YY = Y_coor;
    
    
    for i = 1:size(X_coor,1)
        for j = 1:size(X_coor,2)
            allTest{tag}.X = XX(i,j);
            allTest{tag}.Y = YY(i,j);
            tag = tag+1;
        end
    end
    catch
        1;
    end
end
allTest0 = allTest(1:tag-1);

load('StationaryCheck_UK_bat_387_416.mat')
tag = 1;
for bat = 387:416
    try
    Y_coor = y(3*bat+1 : 3*bat+3);
    X_coor = x(1:end);
    [X_coor,Y_coor] = meshgrid(X_coor,Y_coor);
    XX = X_coor;
    YY = Y_coor;
    
    
    for i = 1:size(X_coor,1)
        for j = 1:size(X_coor,2)
            allTest{tag}.X = XX(i,j);
            allTest{tag}.Y = YY(i,j);
            tag = tag+1;
        end
    end
    catch
        1;
    end
end
allTest0 = [allTest0,allTest];
clear allTest;
%%
[MKpV,SRho,Stat,X,Y] = deal([]);
for i = 1:length(allTest0)
    try
        
        Stat(i) = allTest0{i}.adfH+(1-allTest0{i}.MKendallH)+(1-allTest0{i}.SRhoTd);
        adfpValue(i) = allTest0{i}.adfpValue;
        MKpV(i) = allTest0{i}.MKendallpValue;
        SRho(i) = allTest0{i}.SRhopValue;
        X(i) = allTest0{i}.X;
        Y(i) = allTest0{i}.Y;
    catch
        1;
    end
end
%%

subplot(1,4,1)
data = transpose(buffer(adfpValue,701));
pcolor(data);
set(gca,'YDir','rev');
shading flat;
cptcmap('GMT_drywet', 'mapping', 'scaled');
title('Pval of Dicky-Fuller Test');
colorbar
axis off

subplot(1,4,2)
data = transpose(buffer(SRho,701));
pcolor(data);
set(gca,'YDir','rev');
shading flat;
cptcmap('GMT_haxby', 'mapping', 'direct');
set(gca,'YDir','rev','clim',[0,1]);
title('Pval of Spearman Rho Test');
colorbar
axis off

subplot(1,4,3)
data = transpose(buffer(MKpV,701));
pcolor(data);
set(gca,'YDir','rev');
shading flat;
cptcmap('GMT_haxby', 'mapping', 'direct');
set(gca,'YDir','rev','clim',[0,1]);
title('Pval of Manning Kendall Test')
colorbar
axis off

subplot(1,4,4)
data = transpose(buffer(Stat,701));
pcolor(data);
shading flat;
cbh = colorbar;
cptcmap('GMT_haxby', 'mapping', 'direct');
set(gca,'YDir','rev','clim',[1,3]);
cbh.Ticks = [0,1,2,3] ; %Create 8 ticks from zero to 1
cbh.TickLabels =  {'FullyBad','Bad','Bad','Stationary'};
axis off
title('Stationary Test (1980-2017)')