function [networkProperties,dataOut3,dataOut4,finalStats]=calculateRidgeParams(finalRidges,finalStats,dataIn,mask2D)

%function [networkProperties,dataOut3,dataOut4]=calculateRidgeParams(finalRidges,finalStats,dataIn)
%------- Calculate width of vessels together with other params from scale space analysis of ridges after Lindeberg's algorithm
%------- VARARGIN   :   dataIn              = image to be analysed for ridges on scale space
%-------                finalRidges         = data reduced in dimensions by uniform pyramid 
%-------                finalStats          = [ridgeSaliency; ridgeLength; ridgeWidth; indexSaliency; ridgeWidthCalib]'
%------- Varargout  :   networkProperties   = many parameters: num vessels, tot length, av length, av
%                                               Diameter, av length top 10, relative area ...
%-------                dataOut3            = original image with the ridges overlaid as width 
%-------                dataOut4            = idem but for top 50

%------ no input data is received, error -------------------------
if nargin<1; help calculateRidgeParams;  networkProperties=[];dataOut3=[];dataOut4=[]; return; end;
%%
[rows,cols,levs]                        = size(finalRidges); %#ok<NASGU>


%% Discard those edges covered by a mask
if exist('mask2D','var')
    finalRidges                             = finalRidges.*repmat(mask2D(2:end-1,2:end-1,:),[1 1 levs]);
    fRidges2dB                               = max(finalRidges,[],3);
    fRidges2d                               = zeros(size(fRidges2dB));
    ridgesThatStay                          = unique(finalRidges);
    finalStats                              = finalStats(ridgesThatStay(2:end),:);
    
    [sortedSaliency,indexSaliency]          = (sort(finalStats(:,1)));
    finalStats(:,4)                         = indexSaliency;
    
    %fRidges2d                              = max(finalRidges,[],3);
    %fRidges2d2                             = max(finalRidges2,[],3);
    
    for counterRidges = 1:numel(ridgesThatStay)-1
        fRidges2d                           =fRidges2d+ counterRidges*(fRidges2dB==ridgesThatStay(counterRidges+1));
    
        
        
    end
    %clear finalRidges
%    finalRidges                             = finalRidges2;
 %   clear finalRidges2
else
    fRidges2d                               = max(finalRidges,[],3);
end


%%

if ~exist('pixCalibration','var'); pixCalibration=1.75; end

%% Basic numbers
numRidges                               = size(finalStats,1);
top10                                   = max(1,numRidges-9):numRidges;
top50                                   = max(1,numRidges-49):numRidges-10;

indexTop10                              = finalStats(top10,4);
%indexTop50                              = finalStats(top50,4);

%fRidges2d                               = sum(finalRidges,3);
% To calculate the projection of the ridges in 2D  the  sum can create problems where there are overlaps
% of ridges, it is better to use max and that will assign the overlapped pixel to the deepest ridge.
% the z dim will merge ridges that are not necessarily related.

%% Projection in 2D


%% Correct Length

% to calculate the CORRECT length, get the perimeter of an area, as the area is actually just the
% distance over one of the sides, it does not take into account the slope
actualLength                            = regionprops(fRidges2d,'perimeter','area');

finalStats(:,2)                         = ceil([actualLength(1:numRidges).Perimeter]/2);

%% calculate two separate matrices one with thin&short ridges and one with thick-Long ones.
fRidges2Dthin                           = zeros(rows,cols);
fRidges2Dthic                           = zeros(rows,cols);
fRidges2Dspur                           = zeros(rows,cols);

kThin=1;
kThic=1;


indexThic=[];
indexThin=[];

for k=1:numRidges
    if (finalStats(k,3)>5)&&(finalStats(k,2)>40)
        fRidges2Dthic                   = fRidges2Dthic + (kThic*(fRidges2d==k));
        kThic                           = kThic+1;
        indexThic                       =[indexThic;k]; %#ok<AGROW>
    else
        fRidges2Dthin                   = fRidges2Dthin + (kThin*(fRidges2d==k));
        kThin                           = kThin+1;
        indexThin                       =[indexThin;k]; %#ok<AGROW>
    end
    fRidges2Dspur                       = fRidges2Dspur | bwmorph(fRidges2d==k,'spur',1);
end



%% Correct Width

% For those ridges that are thick (below level 5 of the scale space) and relatively long (more than 30
% pixels long) calculate through the canny edges on the data, otherwise, just adjust from the depth

sizeThinRidges                          = regionprops(fRidges2Dthin,'Area','Centroid','BoundingBox','orientation');
sizeThicRidges                          = regionprops(fRidges2Dthic,'Area','Centroid','BoundingBox','orientation');
%fRidges2DthicLong_L                     = bwlabel(ismember(fRidges2Dthic,find([sizeThicRidges.Area]>40)));
%%
%indexThicLong                           = indexThic([sizeThicRidges.Area]>40);
%%
%sizeThicLongRidges                      = regionprops(fRidges2DthicLong_L,'Area','Centroid','BoundingBox','orientation');

%%
fRidges2d                               = padData(fRidges2d,1);
fRidges2Dthin                           = padData(fRidges2Dthin,1);
fRidges2Dthic                           = padData(fRidges2Dthic,1);
fRidges2Dspur                           = padData(fRidges2Dspur,1);
 %fRidges2DthicLong_L     = padData(fRidges2DthicLong_L,1);


%%
%borders1                = (watershed(bwdist(bwmorph(fRidges2Dthin,'spur',1)|bwmorph(fRidges2Dthic,'spur',1)))==0);
borders1                                = (watershed(bwdist(fRidges2Dspur))==0);

%% Calculate mask for thick/thin ridges plus 

[finalBoundary,fBoundaryStats]          = findVessBoundary(dataIn,fRidges2d,fRidges2Dthic,sizeThicRidges,borders1);
[finalBoundary3]                        = findVessBoundary(dataIn,fRidges2d,fRidges2Dthin,sizeThinRidges,borders1,finalStats(indexThin,:));

if isempty(finalBoundary); finalBoundary=zeros(size(finalBoundary3)); end
if isempty(finalBoundary3); finalBoundary3=zeros(size(finalBoundary)); end

%% recallibrate the values for thick/thin ridges
finalStats(indexThic,6)                 = fBoundaryStats;
finalStats(indexThin,6)                 = 2*finalStats(indexThin,3);

finalStats(finalStats(:,6)==-1,6)       = 2*finalStats(finalStats(:,6)==-1,3);


%% Final Stats from the image
%Number of vessels

if exist('mask2D','var')
    totArea                                     = sum(mask2D(:));
else
    totArea                                     = rows*cols;
end
%%    
numLongVessels                                  = sum(finalStats(:,2)>20);
totLength                                       = sum(finalStats(:,2));

networkProperties.numVessels                    = numRidges;
networkProperties.totLength                     = sum(finalStats(:,2));
networkProperties.avDiameter                    = mean(finalStats(:,6));
networkProperties.avLength                      = mean(finalStats(:,2));

networkProperties.totLength_top10               = sum(finalStats(indexTop10,2));
networkProperties.avDiameter_top10              = mean(finalStats(indexTop10,5));
networkProperties.avLength_top10                = mean(finalStats(indexTop10,2));


networkProperties.numLongVessels                = numLongVessels;
networkProperties.numVesselsPerArea_um2         = numRidges/(rows*pixCalibration)/(cols*pixCalibration);
%networkProperties.numVesselsPerArea_um2         = numRidges/(rows*pixCalibration)/(cols*pixCalibration);
networkProperties.numVesselsPerArea_um2         = numRidges/(totArea*(pixCalibration)^2);
%networkProperties.numLongVesselsPerArea_um2     = numLongVessels/(rows*pixCalibration)/(cols*pixCalibration);
networkProperties.numLongVesselsPerArea_um2     = numLongVessels/(totArea*(pixCalibration)^2);


%Length of vessels
top50length                                     = mean(finalStats(finalStats([top10 top50],4),2));
networkProperties.avLengthTop_um                = pixCalibration*top50length;
%networkProperties.totLengthPerArea_um2          = totLength/rows/cols/pixCalibration;
networkProperties.totLengthPerArea_um2          = totLength/totArea/pixCalibration;



%width of vessels
top50width                                      = mean(finalStats(finalStats([top10 top50],4),6));
networkProperties.avDiameterTop_um              = pixCalibration*top50width;

%Covered area mask for all vessels
%[relAreaCovered,dataOut3]                       = vesselAreaMask(finalRidges,finalStats,dataIn);
%networkProperties.relAreaCovered                = relAreaCovered;

%networkProperties.relAreaCovered                = sum(sum((finalBoundary|finalBoundary3)))/rows/cols;
networkProperties.relAreaCovered                = sum(sum((finalBoundary|finalBoundary3)))/totArea;




%% Calculate the final Masks

[rowsD,colsD,levsD]                             = size(dataIn);

if ~isa(dataIn,'uint8')
    dataIn = uint8(dataIn);
end

dataOut3                                        = dataIn;


dataOut3(:,:,1)                                 = dataIn(:,:,1).*uint8(1-zerocross((finalBoundary|finalBoundary3)-0.5));


% [rowsD,colsD,levsD]                             = size(dataIn);
% dataOut3                                        = dataIn;
% if isa(dataIn,'uint8')
%     dataOut3(:,:,1)                                 = dataIn(:,:,1).*uint8(1-zerocross((finalBoundary|finalBoundary3)-0.5));
% else
%     dataOut3(:,:,1)                                 = dataIn(:,:,1).*(1-zerocross((finalBoundary|finalBoundary3)-0.5));
%     dataOut3                                    =uint8(dataOut3);
%     end


dataOut4                                        = dataIn;
dataOut4(:,:,1)                                 = dataIn(:,:,1).*uint8(1-(finalBoundary|finalBoundary3))+0.8*dataIn(:,:,1).*uint8((finalBoundary|finalBoundary3)).*uint8(1-fRidges2d);
if levsD==1
    dataOut3(:,:,2)                                 = dataIn(:,:,1).*uint8(1-fRidges2d);
    dataOut4(:,:,2)                                 = dataIn(:,:,1).*uint8(1-fRidges2d);
    dataOut3(:,:,3)                                 = dataIn(:,:,1);
    dataOut4(:,:,3)                                 = dataIn(:,:,1);
else
    dataOut3(:,:,2)                                 = dataIn(:,:,2).*uint8(1-fRidges2d);
    dataOut4(:,:,2)                                 = dataIn(:,:,2).*uint8(1-fRidges2d);
end
%%
% if nargout>2
%     %Covered area mask for top 50 vessels
%     if numRidges>50
%         %indexStats                                  = finalStats(max(1,numRidges-49):numRidges,4);
%         %indexStats(indexStats<1)                    = [];
%         %finalRidges                                 = bwlabeln(ismember(finalRidges,[top10 top50])>0);%.*finalRidges;
%         [relAreaCoveredtop50,dataOut4]              = vesselAreaMask(finalRidges,finalStats,dataIn,[indexTop10; indexTop50]);
%         networkProperties.relAreaCoveredtop50       = relAreaCoveredtop50;
%     else
%         dataOut4                                    = dataOut3;
%         networkProperties.relAreaCoveredtop50       = relAreaCovered;
%     end
%     
% end




