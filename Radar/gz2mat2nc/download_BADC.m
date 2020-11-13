clear;
clc;
close all
warning off

addpath(genpath(cd))
mkdir('G:\NIMROD_2019_2020')
cd 'G:\NIMROD_2019_2020'
addpath(genpath(cd))

f = ftp('ftp.ceda.ac.uk','ychen021','AaBb14207');

for year = 2020
    fol = strcat(['/badc/ukmo-nimrod/data/composite/uk-1km/',num2str(year),'/']);
    
    details = dir(f,fol);
    tar_fol = 'G:\NIMROD_2019_2020';
    
    for i = 1:length(details) %%%%%%%%%%%%
        loc = details(i).name;
        mget(f,[fol,loc],tar_fol);
    end
    
end















