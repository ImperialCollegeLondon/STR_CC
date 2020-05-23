
filename = 'C:\Users\Yuting Chen\Downloads\hydrosheds-882eab271075e2e3700e\eu_bas_30s_beta\eu_bas_30s_beta.shp';
S = shaperead(filename);

%%
figure;hold on;

arrayfun(@(s)plot1(s),S,'uniformoutput',false);
hold off
area = arrayfun(@(s)get1(s),S,'uniformoutput',false);

function plot1(s)
if s.BoundingBox(1)<5 & ...
        s.BoundingBox(3) > 49 & s.BoundingBox(4) < 60 & ...
        s.AREA_SQKM > 50 %#ok<AND2>
fill([s.X,s.X(1)],[s.Y,s.Y(1)],'r');
end
end

function a = get1(s)
if s.BoundingBox(1)<5 & ...
        s.BoundingBox(3) > 49 & s.BoundingBox(4) < 60 & ...
        s.AREA_SQKM > 50 %#ok<AND2>
a = s.AREA_SQKM;
else
a = [];
end
end



