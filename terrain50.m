% Read Terrain 50 dtm 
% Source OS survey
% @ Yuting Chen
% yuting.chen17@imperial.ac.uk
% Update: 2020.01.19

clear;clc

% filePath = 'K:\UK_shape\terr50_gagg_gb\data\';
% PATHS = dir([filePath,'*\*_OST50GRID_*.zip']);
% for i = 1:length(PATHS)
%     zipPath = [PATHS(i).folder,filesep,PATHS(i).name];
%     unzip(zipPath, 'K:\UK_shape\terr50_gagg_gb\Terrain50UK')
% end


% filePath = 'K:\UK_shape\terr50_gagg_gb\Terrain50UK\';
% PATHS = dir([filePath,'*.asc']);
% % Preallocate space for DTM50 matrix.
% Eno = 0:50:660000;
% Nno = 1230000:-50:0;
% DTM50 = NaN(length(Nno),length(Eno));
% 
% for i = 1:length(PATHS)
%     
%     FILENAME = [filePath,PATHS(i).name];
%     
%     [OUT, ncols, nrows, xllcorner, yllcorner, cellsize, nodata] = ...
%         ascii_reader(FILENAME,5);
%     if isempty(OUT)
%         [OUT, ncols, nrows, xllcorner, yllcorner, cellsize, nodata] = ...
%         ascii_reader(FILENAME,6);
%     end
%     Eind = findInd(xllcorner,Eno);
%     Nind = findInd(yllcorner+(nrows-1)*cellsize,Nno);
%     DTM50(Nind:Nind+nrows-1,Eind:Eind+ncols-1) = OUT;
%     
% end
% 
% save(['K:\UK_shape\DTM50.mat'],'DTM50','Eno','Nno','-v7.3');

load(['K:\UK_shape\DTM50.mat'],'DTM50','Eno','Nno')

DTM500 = imresize(DTM50,0.1,'bilinear');
Eno = Eno(1):500:Eno(end);
Nno = Nno(1):-500:Nno(end);

pcolor(Eno,Nno,DTM500);
shading flat



function ind = findInd(x,Vect)

ind = find(round(Vect) == x,1);

end



