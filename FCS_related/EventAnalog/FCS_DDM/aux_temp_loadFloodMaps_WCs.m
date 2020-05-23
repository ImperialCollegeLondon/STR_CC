

function [floodmaps,E,N] = aux_temp_loadFloodMaps_WCs(eventNo,unit,mode)
% 
% Example:
%     [floodmaps,E,N] = aux_temp_loadFloodMaps_WCs(3);

arguments
    eventNo (1,:) double
    unit (1,:) char = '1m'
    mode (1,:) char = 'maxFilter'
end

reso = 1;
[ER,NR] = getENRange();
EVec = ER(1):reso:ER(2);
NVec = NR(1):reso:NR(2);


floodmaps = [];
filePath = 'G:\BIGDATA\TOPIC 2\ProcessedFiles\WCMaps_accurate';



for evi = 1:numel(eventNo)
    floodmaps{evi} = sparse([],[],[],numel(EVec),numel(NVec));
    
    for WCNo = 1:19
        try
            fileName = sprintf('WCMaps_Event%03d_WCNo%02d.mat',eventNo(evi),WCNo);
            A = load([filePath,filesep,fileName],'FloodVec','EVec','NVec');
        catch
            continue
        end
        depth = A.FloodVec(A.FloodVec>0);
        A.EVec = A.EVec(A.FloodVec>0);
        A.NVec = A.NVec(A.FloodVec>0);
        if ~isempty(A.EVec)
            [loci,locj] = locateVec(ER,NR,reso,A.EVec,A.NVec);
            % floodmaps(:,:,evi) = A.FMap;
            for pi = 1:numel(depth)
                floodmaps{evi}(loci(pi),locj(pi)) = depth(pi);
            end
        end
    end
end

E = EVec;% A.E;
N = NVec;% A.N;

end


function [ER,NR] = getENRange()
load('G:\BIGDATA\TOPIC 2\ProcessedFiles\FloodMaps_40_30\FloodMaps_Merged_EventNo003.mat','E','N')
ER = [min(E(:)),max(E(:))];
NR = [min(N(:)),max(N(:))];
end

function [loci,locj] = locateVec(ER,NR,reso,E,N)
loci = round((E-ER(1))./reso)+1;
locj = round((N-NR(1))./reso)+1;
end





