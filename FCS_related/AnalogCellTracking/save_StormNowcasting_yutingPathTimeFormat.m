% Extract files for DataDriven Nowcasting
function [] = save_StormNowcasting_yutingPathTimeFormat(cellType)
%
% Saving Format
% Save cell info for each certain time.
% time equal to [1,cell duration]
%
% Data from same interval from cell origin will be stored in the same file
% while the interval from storm origin will be stored as a property for
% each snapshop of each cell path.
%

% cellType = 'SCC';%'LSC';
sfp = 'E:\OneDrive - Imperial College London\RAIN++_OnlyCode&MeetingNotes\PhD_RainPlusPlus\Birmingham_CellPath\';
savefp = 'G:\BIGDATA\TOPIC 2\ProcessedFiles\Nowcasting\';
mkdir(savefp);

STORMFP = getStormFP(sfp,cellType);

try
    
    props = [];
    
    for stormFP = STORMFP'
        tic
        % INITIALIZE
        
        stormFP = stormFP{1};
        stormOrigin = datetime(stormFP(end-14:end-1),'format','yyyyMMddHHmmSS');    [pathNos,pathInfos] = getPathFP(stormFP);
        
        for pathNo = reshape(pathNos,1,[])
            
            thisPath = pathInfos(pathNo,:);
            thisPathProps = pathData(stormFP,cellType,pathNo);
            startT2Storm = minutes(datetime(thisPath.startingTime,'format','dd/MM/yyyy  HH:mm:SS')...
                -stormOrigin)/5;% because radar resolution is 5min;
            
            for interval = 1:thisPath.pathLength
                
                % update for each interval
                
                props_aux = DERIVECELL(thisPathProps,interval,startT2Storm);
                props_aux.isnewstart = thisPath.isNewStart;
                props_aux.traceStormID = str2num(sprintf('%s',stormOrigin));
                props_aux.tracecellTypeID = sprintf('%s',cellType);
                props_aux.tracepathNoID = pathNo;
                
                % props_aux = addTimeInd(props_aux,tstep,TStep{1},cellType);
                try
                    props{interval} = appendT(props{interval},props_aux);
                catch me
                    props{interval} = [];
                    props{interval} = appendT(props{interval},props_aux);
                end
                
            end
            
        end
        
        spp = sprintf('%s%s%s',savefp,filesep,cellType);
        mkdir(spp);
        for tempi = 1:length(props)
            fileName = sprintf('analogFeature_interval%03d',tempi);
            saveName = sprintf('%s%s%s.csv',spp,filesep,fileName);
            if ~isempty(props{tempi})
                thistable = struct2table(props{tempi});
                fid = fopen(saveName);
                if fid~=-1
                    oldTable = readtable(saveName);
                    thistable = [oldTable;thistable]; %#ok<AGROW>
                else
                    %
                end
                writetable(thistable,saveName);
                try
                    fclose(fid);
                catch
                end
            end
        end
        props = []; % clear up and append to exstihng file next time.
        toc
    end
    
catch
    1;
end

end

%% AUXILLARY FUNCTION
function STORMFP = getStormFP(sfp,cellType)
%
% Input: ...
% Output: STORMFP: <cell> eg.'...\2005\LSC\20050521133500\'
%

stormDir = dir([sfp,'*\',cellType,'\201*']);
STORMFP = arrayfun(@(x)[x.folder,'\',x.name,'\'],stormDir,...
    'UniformOutput',false);

end

function [pathNos,pathInfos] = getPathFP(stormFP)
%
% all cell paths are taken into consideration
% Output: PathNos: [:,1]
pathSummary = [stormFP,'cellPathSummary.csv'];
T = readtable(pathSummary);
pathNos = 1:size(T,1);
pathInfos = T;
end

function CIND = getCellID(stormFP,cellType,tstep)
%
% Output: CIND: <int> eg. 81
%
fileN = [stormFP,'\Property\',tstep,'_Property_',cellType,'.csv'];
CIND = size(readtable(fileN),1);

end

function props_aux = addTimeInd(props_aux,tstep,startT,cellType)
%
% Output: additional columes in table for labeling time.
%

props_aux.tstep = repmat(tstep,size(props_aux,1),1);

ctid = [1:size(props_aux,1)]';
props_aux.TRACEID = func111(ctid,tstep,startT,cellType);

    function traceid = func111(ctid,tstep,startT,cellType)
        
        traceid = {};
        for i = ctid'
            traceid{i,1} = {sprintf('%s_%s_%s_%d',startT,cellType,tstep,ctid(i))};
        end
        
    end

end

function stormT = appendT(stormT,props_aux)
% ......
% think about how to make it faster.
% ......

stormT = [stormT;props_aux];

end

function thisPathProps = pathData(stormFP,cellType,pathNo)

% GET ALL RELATED INFO FOR THIS PATH

thisPathProps = struct;

thisPathProps.Ref = readtable(sprintf('%scellPath_Ref%sRef_path%d.csv',...
    stormFP,filesep,pathNo),'Delimiter',',');
RCell = String2DoubleV(thisPathProps.Ref.Reflectivity);
thisPathProps.Ref.Intensity = getInt_in_RCell(RCell);
thisPathProps.Prop = readtable(sprintf('%scellPath_Prop%sProp_path%d.csv',...
    stormFP,filesep,pathNo));
thisPathProps.Pixels = readtable(sprintf('%scellPath_Pixels%sPixels_path%d.csv',...
    stormFP,filesep,pathNo),'Delimiter',',');

end



function props_int = DERIVECELL(thisPathProps,interval,startT2Storm)
%
% DERIVE characteristics for an interval for pathNo in stormFP celltype.

% PropertyName:

try
    props_aux = thisPathProps.Prop;
    props_aux.IMF = cell2mat(cellfun(@(x)nanmean(x),thisPathProps.Ref.Intensity,...
        'UniformOutput',false));
    % prope_aux = join(props_aux,..)
    
    % calculate props_int for this interval
    
    IND = 1:interval;
    props_int = struct;
    props_int.totalDepth = nansum(props_aux.IMF(IND));
    props_int.medianArea = nanmedian(props_aux.Area(IND));
    props_int.velmag = nanmedian(props_aux.Velocity(IND));
    props_int.velori = nanmedian(props_aux.Orientation(IND));
    props_int.skewness = skewness(props_aux.IMF(IND));
    props_int.currentInt = props_aux.IMF(IND(end));
    props_int.dur2storm = interval+startT2Storm-1;
    
    
    % climate related
    % need to be updated
    props_int.humid = 0;
    props_int.temper = 0;
    props_int.windS = 0;
    props_int.seaLP = 0;
    
    
    %{
    [Id,Area,Centroid_1,Centroid_2,WeightedCentroid_1,WeightedCentroid_2,Perimeter,...
        AxisMaj,AxisMin,Eccen,Orientation,NumCell,MImaj,MImin] = importProps_1(stormFP,cellType,tstep,cInd);
    
    [Velocity,Angle] = importProps_2(stormFP,cellType,tstep,cInd);
    
    [Rmean, Rmax, Istd, IMin, IMean, IMax, StormX, StormY] = importProps_3(stormFP,cellType,tstep,cInd);
    
    props = table(Id,WeightedCentroid_1,WeightedCentroid_2,Centroid_1,Centroid_2,Area,Perimeter,...
        IMax, IMin, IMean,...
        Orientation,AxisMaj,AxisMin,Eccen,...
        Rmean, Rmax, Velocity,Angle,NumCell,...
        Istd, MImaj,MImin,StormX, StormY);
    %}
catch
    1;
end

    function [Id,Area,Centroid_1,Centroid_2,WeightedCentroid_1,WeightedCentroid_2,Perimeter,AxisMaj,AxisMin,Eccen,Orientation,NumCell,MImaj,MImin] = importProps_1(stormFP,cellType,tstep,cInd)
        % import Props from '/Property'
        % 1.Area;
        % 3.Centroid_1; Centroid_2;
        % 6.AxisMaj; AxisMin; Orientation;
        % 7.NumCell;
        % 5.MImaj; MImin;
        fileN = [stormFP,'\Property\',tstep,'_Property_',cellType,'.csv'];
        T = readtable(fileN);
        Area = T.Area; % Threshold: 35 dBz for LSC
        %
        % #NOTICE# make sure Centroid is in reference to the 'left right corner'
        %
        Id = T.Id;
        Centroid_1 = T.Centroid_1;
        Centroid_2 = T.Centroid_2;
        WeightedCentroid_1 = T.WeightedCentroid_1;
        WeightedCentroid_2 = T.WeightedCentroid_2;
        Perimeter = T.Perimeter;
        AxisMaj = T.AxisMaj;
        AxisMin = T.AxisMin;
        Eccen = T.Eccen;
        Orientation = T.Orientation;
        %
        % #NOTICE# make sure which the reference of the orientation is
        %
        NumCell = repmat(cInd,size(AxisMaj));
        [MImaj,MImin] = computeMI(AxisMaj,AxisMin);
        
        function [MImaj,MImin] = computeMI(AxisMaj,AxisMin)
            % 5. Moment of Inertia in both major Axis and minor Axis
            % Reference: https://en.wikipedia.org/wiki/List_of_second_moments_of_area
            a = AxisMaj/2;
            b = AxisMin/2;
            MImaj = pi/4*a.*b.^3;
            MImin = pi/4*a.^3.*b;
        end
    end

    function [Velocity,Angle] = importProps_2(stormFP,cellType,tstep,cInd)
        % import Props from '/Velocity'
        % 4.Velocity; Angle;
        fileN = [stormFP,'\Velocity\',tstep,'_Velocity_',cellType,'.csv'];
        T = readtable(fileN);
        Velocity = T.Velocity;
        Angle = T.Angle;
    end

    function [Rainmean, Rainmax, Istd, dbzmin,dbzmean, dbzmax, StormX, StormY] = importProps_3(stormFP,cellType,tstep,cInd)
        % import Props from '/Property' + '/Reflectivity'
        % 2.Imean; Imax; dbzmean; Istd; dbzmax;
        % 8.StormX; StormY
        
        % from '/Property'
        % Imean; Imax; dbzmean; dbzmax;
        fileN = [stormFP,'\Property\',tstep,'_Property_',cellType,'.csv'];
        T = readtable(fileN);
        dbzmin = T.IMin;
        dbzmean = T.IMean;
        dbzmax = T.IMax;
        Rainmean = getInt(dbzmean);
        Rainmax = getInt(dbzmax);
        
        % from '/Reflectivity'
        % Istd; StormX; StormY
        fileN = [stormFP,'\Reflectivity\',tstep,'_Reflectivity_',cellType,'.csv'];
        T = readtable(fileN,'Delimiter',',');
        RCell = String2DoubleV(T.Reflectivity);
        ICell = getInt_in_RCell(RCell);
        Istd = cellfun(@(x)std(x),ICell);
        
        % .....
        StormX = repmat(500,size(Istd)); StormY = repmat(500,size(Istd));
        %
        % #NOTICE# here the stormX/Y need further discussion.
        %
        % .....
        
        
    end

end

function RCell = String2DoubleV(StringT)
% TRANSFORM from cells of strings to cell of double vector.
% Output: RCell n*1 cell. each cell is a double vector (1*m)
RCell = cellfun(@(x)str2num(x),StringT,'UniformOutput',false); %#ok<ST2NM>
end

function Imean = getInt(dbzMean)
% dBz 2 intensity;
% format: all saved in double vector.
a = 200;
b = 1.6;
Imean = (10.^(dbzMean/10)/a).^(1/b);
end

function ICell = getInt_in_RCell(RCell)
% dBz 2 intensity;
% format: all saved in cell.
%
a = 200;
b = 1.6;
ICell = cellfun(@(x)(10.^(x/10)/a).^(1/b),RCell,'UniformOutput',false);

end



% 1. Area - could you also specify which threshold you used for recognition?
% LSC: 35dbz
% Resolution: 1 km * 1 km
% // get
% // from '/Property' .Area
%
% 3. Centroid Position - in reference to the left right corner of the study area mask would be perfect, but absolute values can also work.
% // get
% centroid position is in reference to the left right corner: fine
% // from  '/Property' .Centroid_1, .Centroid_2  (WeightedCentroid_1/2 might can be used as well)
%
% 6. Ellipse Characteristics : Length of Major & Minor Axis, angle of the Major Axis
% // get
% // from 'Property' .AxisMaj, .AxisMin, .Orientation
%
% 7. Number of Storm Cells - so i can account for the splitting and the merging of the timesteps
% Means the total cell number in this time step
% // count total cells in this time step
% // from 'Property' num($Id$)
%
% 5. Moment of Inertia in both major Axis and minor Axis
% ? choose a tool to derive
% ? x y direction means x/y of elispse
% https://archive.cnx.org/contents/2cd75bcc-817d-44c7-9d8b-cec23856d26d@1/implementation-of-moment-of-inertia
% pixels
% // from 'Property' as well.
%
%
% 4. Velocity & Direction  of the storm centroid. You can either use velocity in m/s and direction or pixel displacement for each X and Y direction. Either way could you also specify the +/- senses?
% It is velocity and angle in property file
% // get
% // from '/Velocity' .Velocity, .Angle
%
%
% 2. Intensity : mean, maximum, and standard deviation values
% Statistics in space.
% // get
% // mean,maximum: from '/Property'
% // standard deviation values: from '/Reflectivity' then Compute
%
%
% 8. Storm Extent in X and Y Direction
% They mean area or what extent in X and Y direction ????? %%% waiting for reply from Lipen
% Regarding the storm extent, it is simply the extension of the storm (length in km or number of radar cells) in the x and y direction.
% Xmin, Ymin, Xmax and Ymax in other words.
% // still puzzeled.
% // from '/Reflectivity'
