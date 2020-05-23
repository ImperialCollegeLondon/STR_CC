clear;
clc;
close all
warning off

addpath(genpath(cd))
mkdir('F:\NIMROD_2016_2019')
cd 'F:\NIMROD_2016_2019'
addpath(genpath(cd))

f = ftp('ftp.ceda.ac.uk','ychen021','AaBb14207');

for year = 2018
    fol = strcat(['/badc/ukmo-nimrod/data/composite/uk-1km/',num2str(year),'/']);
    
    details = dir(f,fol);
    tar_fol = 'F:\NIMROD_2016_2019';
    
    for i = 221:length(details)
        loc = details(i).name;
        mget(f,[fol,loc],tar_fol);
    end
    
end















