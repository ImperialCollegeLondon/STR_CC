% Read RPS Flood Map
% Notes
% It takes a quite while (10 min! how) just to open one Ascii file (~200Mb)
% 52.4 GB * 4 in total
% @ Yuting Chen
% Update: 2020.01.08

filename = ['G:\BIGDATA\TOPIC 2\FileFromLIPEN\Flood Maps\Batch_1\',...
    '035_WC13_WDraster_PEAK.asc'];    %'All_Batch_1_Output_Merged.asc'];   
% filename = 'G:\BIGDATA\TOPIC 2\FileFromLIPEN\Flood Maps\Batch_1\036_WC14_WDraster_PEAK.asc';
% filename = 'G:\BIGDATA\TOPIC 2\FileFromLIPEN\Flood Maps\Batch_1\115_WC19_WDraster_PEAK.asc';
[OUT, ncols, nrows, xllcorner, yllcorner, cellsize, nodata] = ...
    ascii_reader(filename);

%%
OUT1 = exp(OUT)-1;
imagesc(imresize(OUT,1/250));
box on
%%
out = unique(diff(sort(OUT)));

out2 = sparse(OUT);

S2 = getfield(whos('out2'),'bytes');

S = getfield(whos('out'),'bytes');


