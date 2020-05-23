clear;clc

addpath('C:\Users\Yuting Chen\Dropbox (Personal)\MatLab_Func');

for YEAR = 2018%2007:2016

    [status] = aggregate_NIMROD_Hour(YEAR);
    status

end



