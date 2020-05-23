% ------------------------------------------------------- %
%
% Numerical Experiment for optimizing model configuration
%
% @ yuting
% yuting.chen17@imperial.ac.uk
%
% ------------------------------------------------------- %

%% C1U5-1
clear;clc
[version,testConfig,filefolder] = getExpNoInfo('C1U5-1');
STATS = Main_OfflineCrossValidation(version,testConfig,filefolder);

%% C1U10-1
clear;clc
[version,testConfig,filefolder] = getExpNoInfo('C1U10-1');
STATS = Main_OfflineCrossValidation(version,testConfig,filefolder);
%% C1P5-1
[version,testConfig,filefolder] = getExpNoInfo('C1P5-1');
STATS = Main_OfflineCrossValidation(version,testConfig,filefolder);

%% C1P5-Rad1
clear;clc
[version,testConfig,filefolder] = getExpNoInfo('C1P5-Rad1');
STATS = Main_OfflineCrossValidation(version,testConfig,filefolder);

%% C1P5In-Rad1
clear;clc
[version,testConfig,filefolder] = getExpNoInfo('C1P5In-Rad1');
STATS = Main_OfflineCrossValidation(version,testConfig,filefolder);

%% C1P5-2
[version,testConfig,filefolder] = getExpNoInfo('C1P5-2');
STATS = Main_OfflineCrossValidation(version,testConfig,filefolder);

%% C1P5In-1
clear;clc
[version,testConfig,filefolder] = getExpNoInfo('C1P5In-1');
STATS = Main_OfflineCrossValidation(version,testConfig,filefolder);

%% C1U5-Rad1
clear;clc
[version,testConfig,filefolder] = getExpNoInfo('C1U5-Rad1');
STATS = Main_OfflineCrossValidation(version,testConfig,filefolder);

% C1U10-Rad1
clear;clc
[version,testConfig,filefolder] = getExpNoInfo('C1U10-Rad1');
STATS = Main_OfflineCrossValidation(version,testConfig,filefolder);



%% Also consider: 
% Estimate the time when the rainfall event reaching the peak
% use it to determine the corresponding prioTime we will use for allocating climate data.


%% C2P5-Rad1
clear;clc
[version,testConfig,filefolder] = getExpNoInfo('C2P5-Rad1');
STATS = Main_OfflineCrossValidation(version,testConfig,filefolder);






