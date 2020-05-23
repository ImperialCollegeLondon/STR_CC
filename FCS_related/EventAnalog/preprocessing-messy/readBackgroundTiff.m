
% read one openmap shapefile
% Birm


% tiff format
figure;
set(gcf, 'Position', get(0, 'Screensize'));
readOneTile('sp','06')
title('sp-10m')
% readOneTile('sp','08')
% readOneTile('sp','07')

% minisc_gb format
% readOneMinisc()



function readOneMinisc()
filePath = 'K:\UK_shape\OS_MiniScale_GB\minisc_gb\data\RGB_TIF_compressed\';
FILES = dir([filePath,'MiniScale_(relief2)_R22.tif']);

t = Tiff([FILES(1).folder,'\',FILES(1).name],'r');

imageData = read(t);
imshow(imageData) %'Colormap', summer(round((rand*200))));
hold on;
drawnow

end



function readOneTile(area,num)

filePath = 'K:\UK_shape\OC_UKOpenMap - Local\';
FILES = dir([filePath,'omlras_gtfc_',area,'\Data\',upper(area),num,'*.tif']);

for fi = 1:numel(FILES)
    
    %%%% need to find out loaction of the tile.
    
    t = Tiff([FILES(fi).folder,'\',FILES(fi).name],'r');
    strc = imfinfo([FILES(fi).folder,'\',FILES(fi).name]);
    x = strc.ModelTiepointTag(4);
    y = strc.ModelTiepointTag(5);
    imageData = read(t);
    unit = 10;%50;
    imageData = imresize(imageData,'scale',1/unit,'method','box');

    x = x:unit:x+unit*(size(imageData,1)-1);
    y = y:-unit:y-unit*(size(imageData,1)-1);
    imshow(imageData,'XData',x,'YData',y) %'Colormap', summer(round((rand*200))));
%     set(gcf, 'Position', get(0, 'Screensize'));
    hold on;
    % axis equal
    colormap(gca, gray(256))
    drawnow
end
axis normal

end




