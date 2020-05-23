function Area = getBoundaryShp(regionName)

switch(regionName)
    case 'Birm'
        fileName = 'K:\UK_shape\bdline_essh_gb\Data\GB\district_borough_unitary_region.shp';
        area = shaperead(fileName);%,'Attributes',{'BoundingBox','X','Y','Name'});
        Area = area(109); % name: {'Birmingham District (B)'}
end
end

%%
% fileName = 'K:\UK_shape\bdline_essh_gb\Data\GB\district_borough_unitary_region.shp';
% area = shaperead(fileName);
% Area = area(109); % name: {'Birmingham District (B)'}
% plot(Area.X,Area.Y,'-');
% 
% expNo = 'C1P5In-Rad1';eventNo4val = 90;evtPlot = 70;
% [version,testConfig,filefolder] = getExpNoInfo(expNo);
% load([filefolder,sprintf('STATS_%s_ev%03d',version,eventNo4val)],'FOREOUTPUT');
% 
% E = FOREOUTPUT.E(:);
% N = FOREOUTPUT.N(:);
% 
% in = inpolygon(E(:),N(:),Area.X,Area.Y);
% 
% figure;
% hold on;
% plot(Area.X,Area.Y,'r-');
% plot(E(in),N(in),'rx');
% plot(E,N,'k.');
% hold off;
% 
% Area.X = E(in);
% Area.Y = N(in);
% shapewrite(Area,'Birmingham.shp');
% 
% A = shaperead('C:\Users\Yuting Chen\Desktop\Birmingham_shapefile\Birmingham.shp');
% plot(A.X,A.Y,'gx');






