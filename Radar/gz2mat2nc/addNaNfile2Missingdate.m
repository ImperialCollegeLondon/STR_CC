year = 2018;
month = 5;
day = 10;
savePath = 'K:\UK_Radar_Matlab';
sf = sprintf('%s%s%04d_%d_%d.mat',savePath,filesep,year,mon,day);

DAT = struct;
for hour = 0:1:23
    for minu = 0:5:55
        timeStr = sprintf('%04d%02d%02d%02d%02d',year,mon,day,hour,minu);
        thisM = struct;
        [thisM.int_gen_hd, thisM.rl_gen_hd, thisM.rl_datsp_hd, char_hd, ...
            thisM.int_datsp_hd, thisM.rr] = deal(NaN);
        eval(['DAT.d',timeStr,'=thisM;']);
    end
end
save(sf,'DAT','-v7.3');