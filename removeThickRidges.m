function [finalRidges,finalStats,networkProperties]=removeThickRidges(fRidges,fStats,netP) %#ok<INUSD>


ridgeLength                             = fStats(:,2);
ridgeWidthCalib                         = fStats(:,5)/2;
ridgeSaliency                           = fStats(:,1);
ridgeWidth                              = fStats(:,3);
%before determining the order of the saliency, Discard ALL ridges that are too short but very wide
ridgesNotTooWide                        = (ridgeLength./(ridgeWidthCalib))>1.2;

[rows,cols,levs]                        = size(fRidges); %#ok<NASGU>

%
[sortedSaliency,indexSaliency]          = (sort(ridgeSaliency(ridgesNotTooWide)));
%%
finalStats                              = [ridgeSaliency(ridgesNotTooWide) ridgeLength(ridgesNotTooWide) ridgeWidth(ridgesNotTooWide) indexSaliency 2*ridgeWidthCalib(ridgesNotTooWide)];
finalRidges                             = bwlabeln(ismember(fRidges,find((ridgesNotTooWide))));



%%
if exist('netP','var')
    if ~exist('pixCalibration','var'); pixCalibration=1.75; end
    %%
    numRidges                               = size(finalStats,1);
    top10                                   = max(1,numRidges-9):numRidges;
    top50                                   = max(1,numRidges-49):numRidges-10;


       
    %% Final Stats from the image
    %Number of vessels
    numLongVessels                                  = sum(finalStats(:,2)>20);
    totLength                                       = sum(finalStats(:,2));

    networkProperties.numVessels                    = numRidges;
    networkProperties.totLength                     = sum(finalStats(:,2));
    networkProperties.avDiameter                    = mean(finalStats(:,5));
    networkProperties.avLength                      = mean(finalStats(:,2));

    networkProperties.totLength_top10               = sum(finalStats(finalStats(top10,4),2));
    networkProperties.avDiameter_top10              = mean(finalStats(finalStats(top10,4),5));
    networkProperties.avLength_top10                = mean(finalStats(finalStats(top10,4),2));


    networkProperties.numLongVessels                = numLongVessels;
    networkProperties.numVesselsPerArea_um2         = numRidges/(rows*pixCalibration)/(cols*pixCalibration);
    networkProperties.numLongVesselsPerArea_um2     = numLongVessels/(rows*pixCalibration)/(cols*pixCalibration);


    %Length of vessels
    top50length                                     = mean(finalStats(finalStats([top10 top50],4),2));
    networkProperties.avLengthTop_um                = pixCalibration*top50length;
    networkProperties.totLengthPerArea_um2          = totLength/rows/cols/pixCalibration;

    %width of vessels
    top50width                                      = mean(finalStats(finalStats([top10 top50],4),5));
    networkProperties.avDiameterTop_um              = pixCalibration*top50width;

    [relAreaCovered]                                = vesselAreaMask(finalRidges);
    networkProperties.relAreaCovered                = relAreaCovered;
end