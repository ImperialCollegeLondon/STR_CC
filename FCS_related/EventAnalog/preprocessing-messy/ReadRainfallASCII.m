% Read Rainfall ASCII file
% Notes
% @ Yuting Chen
% Update: 2020.01.08

filename = ['G:\BIGDATA\TOPIC 2\FileFromLIPEN\Rainfall_in_157Floods\KED_ASCII\',...
    '201209110010.ASC'];    
[OUT, ncols, nrows, xllcorner, yllcorner, cellsize, nodata] = ...
    ascii_reader(filename);
