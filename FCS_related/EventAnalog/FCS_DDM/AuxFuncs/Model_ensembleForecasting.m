function [FOR_test,FOR_analog] = Model_ensembleForecasting(RES_analog,AnalogInput,LIB)


RES_test = AnalogInput.RIMAGESval;
RES_testtime = AnalogInput.raintimeval;
LIB_time = AnalogInput.raintime;
LIB_images = AnalogInput.RIMAGES;


[FOR_test,FOR_analog] = deal(struct('eventNo',[],'time_start',[],'fore_Images',[],...
    'unit','1km-5min','forehour',1,'leadhour',0));
FOR_test.fore_Images = NaN(0,34,30);%[idx,[Loc]]
FOR_analog.fore_Images = NaN(0,34,30,12);%[idx,[loc],ensemble];

for testInd = 1:size(RES_test,3)
    
    idx = testInd;

    forehour = FOR_test.forehour;
    leadhour = FOR_test.leadhour;
    fore_testtime0 = RES_testtime(idx)+(leadhour)/24;
    fore_testtime1 = RES_testtime(idx)+(leadhour+forehour)/24;
    fore_libtime0 = LIB_time(RES_analog(idx,:))+(leadhour)/24;
    fore_libtime1 = LIB_time(RES_analog(idx,:))+(leadhour+forehour)/24;
    fore_anaimages = NaN(34,30,size(RES_analog,2));
    
    for ensNo = 1:size(RES_analog,2)
        
        fore_anaimages(:,:,ensNo) = nansum(LIB_images(:,:,find(LIB_time == fore_libtime0(ensNo)):find(LIB_time == fore_libtime1(ensNo))),3);
        
    end
    fore_testimages = nansum(RES_test(:,:,find(RES_testtime == fore_testtime0):find(RES_testtime == fore_testtime1)),3);
    

    % Attribute to FOR_test/analog struct
    FOR_test.time_start(testInd) = datenum(fore_testtime0);
    FOR_analog.time_start(testInd,:) = datenum(fore_libtime0');
    
    FOR_test.fore_Images(testInd,:,:) = fore_testimages;
    FOR_analog.fore_Images(testInd,:,:,:) = fore_anaimages;
     
end

FOR_test.eventNo = LIB.FloodLib.test.eventNo;
FOR_analog.eventNo = LIB.FloodLib.train.eventNo(RES_analog);

end
