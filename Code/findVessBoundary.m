function [finalBoundary,fBoundaryStats]=findVessBoundary(dataIn,fRidges2D,fRidges2DthicLong_L,sizeThicLongRidges,borders1,fStats)


% function [finalBoundary]=findVessBoundary(dataIn,fRidges,fStats)
%------- Calculate width of vessels by finding Canny's edges in the vicinity of the ridge. It
%------- is only used for a reduced number of ridges (Thick ones) and called from calculateRidgeParams
%------- VARARGIN   :   dataIn              = image to be analysed for ridges on scale space
%-------                fRidges             = data in 2D with labelled ridges 
%-------                fStats              = [ridgeSaliency; ridgeLength; ridgeWidth; indexSaliency; ridgeWidthCalib]'
%------- Varargout  :   finalBoundary       = The boundary for all the ridges as an area idem but for top 50
%-------                fBoundaryStats      = the av width of each ridge analysed

%------ no input data is received, error -------------------------
if nargin<1; help findVessBoundary;  finalBoundary=[];fBoundaryStats=[]; return; end;

[rows,cols,levs]        =  size(dataIn); %#ok<NASGU>

%%

sizeF                   = 9;
filtg                   = gaussF(sizeF,sizeF,1);
dataInSm                = double(imfilter(dataIn,filtg,'replicate'));

finalBoundary           = zeros(rows,cols);

numThicLongR            = numel(sizeThicLongRidges);%#ok<NASGU>


if numThicLongR==0
    finalBoundary=[];fBoundaryStats=[];
else
    %%

    dilFactor3{60}=0;
    for k=60:-1:1;
        dilFactor3{k}   =strel('disk',k,0);
        %dilFactorSize(k)=1+k*2;%size(getnhood(dilFactor3{k}),1);
    end


    if ~exist('fStats','var')
        %fStats was not passed as a variable, the widths and the mask are calculated for each ridge by
        %finding Canny's edges
        fBoundaryStats(numThicLongR)= 0;

        %currRidge=30;
        dilFactor0  = 45;
        for currRidge=1:numThicLongR
            currRowsD = max(1,ceil(sizeThicLongRidges(currRidge).BoundingBox(2))-dilFactor0):min(rows,floor(sizeThicLongRidges(currRidge).BoundingBox(2))+(sizeThicLongRidges(currRidge).BoundingBox(4))+dilFactor0);
            currColsD = max(1,ceil(sizeThicLongRidges(currRidge).BoundingBox(1))-dilFactor0):min(cols,floor(sizeThicLongRidges(currRidge).BoundingBox(1))+(sizeThicLongRidges(currRidge).BoundingBox(3))+dilFactor0);

            if currRidge==58
                qqq=1;
            end
            %finalBoundary(currRowsD,currColsD)= finalBoundary(currRowsD,currColsD) | findVessBoundary2(dataInSm(currRowsD,currColsD),fRidges2D(currRowsD,currColsD),fRidges2DthicLong_L(currRowsD,currColsD),sizeThicLongRidges,27,rows,cols);
            try
                [tempBoundary,fBoundaryVal]         = findVessBoundary2(dilFactor3,dataInSm(currRowsD,currColsD),fRidges2D(currRowsD,currColsD),fRidges2DthicLong_L(currRowsD,currColsD),sizeThicLongRidges,currRidge,borders1(currRowsD,currColsD));
                finalBoundary(currRowsD,currColsD)  = finalBoundary(currRowsD,currColsD) | tempBoundary;
                fBoundaryStats (currRidge)          = fBoundaryVal ;

            catch
                fBoundaryStats (currRidge)          = -1 ;
                qqq=1; %#ok<NASGU>
            end
        end
    else
        %fStats exist, these are the small ridges that will be uniformly dilated
        dilFactor0  = 45;
        for currRidge=1:numThicLongR
            currRowsD = max(1,ceil(sizeThicLongRidges(currRidge).BoundingBox(2))-dilFactor0):min(rows,floor(sizeThicLongRidges(currRidge).BoundingBox(2))+(sizeThicLongRidges(currRidge).BoundingBox(4))+dilFactor0);
            currColsD = max(1,ceil(sizeThicLongRidges(currRidge).BoundingBox(1))-dilFactor0):min(cols,floor(sizeThicLongRidges(currRidge).BoundingBox(1))+(sizeThicLongRidges(currRidge).BoundingBox(3))+dilFactor0);

            %finalBoundary(currRowsD,currColsD)= finalBoundary(currRowsD,currColsD) | findVessBoundary2(dataInSm(currRowsD,currColsD),fRidges2D(currRowsD,currColsD),fRidges2DthicLong_L(currRowsD,currColsD),sizeThicLongRidges,27,rows,cols);
            [tempBoundary]                      = findVessBoundary3(dilFactor3,fRidges2DthicLong_L(currRowsD,currColsD)==currRidge,fStats(currRidge,3));
            finalBoundary(currRowsD,currColsD)  = finalBoundary(currRowsD,currColsD) | tempBoundary;

            %         figure(5)
            %         surfdat(finalBoundary)
            %         drawnow
            %         pause(0.05)
        end


    end
end

%numRidges               = size(fStats,1);

%fRidges2DthicLong_L       = bwlabel(ismember(fRidges2Dthic,find([sizeThicRidges.Area]>40)));
%sizeThicLongRidges      = regionprops(fRidges2DthicLong_L,'Area','Centroid','BoundingBox','orientation');

%numThicR                = numel(sizeThicRidges); 
%numThinR                = numel(sizeThinRidges);



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [finalBoundary,fBoundaryVal]   = findVessBoundary2(dilFactor3,dataInSm,fRidges2D,fRidges2DthicLong_L,sizeThicLongRidges,currRidge,borders1)
%This function is for those ridges that are THICK and therefore a special calculation of width will be
%performed, Canny edge detection is applied to the vicinity of the  ridge and distance from the ridge to
%the edges on both sides of the ridge is calculated

%%
dilFactor0  = 45;
dilFactor   =strel('disk',dilFactor0,0);
dilFactor2  =strel('disk',5,0);

[rows,cols] =   size(dataInSm);



%%

%currRowsD = max(1,ceil(sizeThicLongRidges(currRidge).BoundingBox(2))-dilFactor0):min(rows,floor(sizeThicLongRidges(currRidge).BoundingBox(2))+(sizeThicLongRidges(currRidge).BoundingBox(4))+dilFactor0);
%currColsD = max(1,ceil(sizeThicLongRidges(currRidge).BoundingBox(1))-dilFactor0):min(cols,floor(sizeThicLongRidges(currRidge).BoundingBox(1))+(sizeThicLongRidges(currRidge).BoundingBox(3))+dilFactor0);
%%
currRidgeIm             = (fRidges2DthicLong_L==currRidge);
currData                = dataInSm.*imdilate(currRidgeIm,dilFactor);
currData2               = currData;
currData2(currData==0)  = 999;
minData2                = min(currData2(:))-1;
currData2(currData2==999)= minData2;
%%
%currRows = (ceil(sizeThicLongRidges(currRidge).BoundingBox(2))):(floor(sizeThicLongRidges(currRidge).BoundingBox(2))+(sizeThicLongRidges(currRidge).BoundingBox(4)));
%currCols = (ceil(sizeThicLongRidges(currRidge).BoundingBox(1))):(floor(sizeThicLongRidges(currRidge).BoundingBox(1))+(sizeThicLongRidges(currRidge).BoundingBox(3)));

%% Find End Points of the Ridge, 
% Calculate the orientation of the last section for:
%    (a) dividing the edge into 2 regions for calculating the asymetric width and
%    (b) calculating a stop line to prevent thick ridges having a round end (THIS IS USED BELOW)

clear angle1 angle2;
%%

numPointsRidge                   = ceil(sum(currRidgeIm(:))/3);
edgesCurrRidge                  = EndPoints(bwmorph(currRidgeIm,'spur',min(6,numPointsRidge)));
%locate the position
indexCurrEdges                  = find (edgesCurrRidge);
%create two matrices, each with one of the end points
borders41                       = zeros(size(currRidgeIm));
borders42                       = zeros(size(currRidgeIm));
borders41(indexCurrEdges(1))    = 1;
borders42(indexCurrEdges(2))    = 1;
%%
%trim the ridge to the last 5-10 points to get the orientation of the ridge CLOSE to the endpoint
bordersEdge1                    = bwmorph(currRidgeIm,'spur',min(6,numPointsRidge)).*imdilate(borders41,dilFactor3{5});
bordersEdge2                    = bwmorph(currRidgeIm,'spur',min(6,numPointsRidge)).*imdilate(borders42,dilFactor3{5});

%%
angle1                          = regionprops(bwlabel(bordersEdge1),'orientation');
angle2                          = regionprops(bwlabel(bordersEdge2),'orientation');

% for the extension line, do not add anything
EndPoint1_extension             = strel('line', 280,angle1(1).Orientation );
EndPoint2_extension             = strel('line', 280,angle2(1).Orientation);

currRidgeImDil                  = imdilate(currRidgeIm,ones(5));
extendedEdge1                  = imdilate(borders41,EndPoint1_extension);
extendedEdge2                  = imdilate(borders42,EndPoint2_extension);

%
[extendedEdge11,numExEdges11]  = bwlabel((extendedEdge1.*(1-currRidgeImDil)));%.*(1-imdilate(extendedEdge2,[1 1;1 1]))));
[extendedEdge22,numExEdges22]  = bwlabel((extendedEdge2.*(1-currRidgeImDil)));%.*(1-imdilate(extendedEdge1,[1 1;1 1]))));

% figure(1)
% subplot(131);surfdat(extendedEdge1+4*currRidgeIm)
% subplot(132);surfdat(extendedEdge2+4*currRidgeIm)
% subplot(133);surfdat(extendedEdge22+2*extendedEdge11+4*currRidgeIm)

%% Remove fragmented edges that DO NOT touch the edges of the image

for counterEdge =1:numExEdges11
    [rrr,ccc]= find(extendedEdge11==counterEdge);
    %test rows first then columns
    if ~any([ismember([1 rows],rrr) ismember([1 cols],ccc)])
        extendedEdge11(extendedEdge11==counterEdge)=0;
    end
end
for counterEdge =1:numExEdges22
    [rrr,ccc]= find(extendedEdge22==counterEdge);
    %test rows first then columns
    if ~any([ismember([1 rows],rrr) ismember([1 cols],ccc)])
        extendedEdge22(extendedEdge22==counterEdge)=0;
    end
end

% 
% figure(2)
% subplot(131);surfdat(extendedEdge11+4*currRidgeIm)
% subplot(132);surfdat(extendedEdge22+4*currRidgeIm)
% subplot(133);surfdat(extendedEdge22+2*extendedEdge11+4*currRidgeIm)
%% now find the two extended edges closest to the ridge

[extendedEdge11,numExEdges11] = bwlabel(extendedEdge11);
[extendedEdge22,numExEdges22] = bwlabel(extendedEdge22);

%
%qqq=bwdist(currRidgeIm);

%% careful examination found that the order in which the edges are found can change due to the spur...

edgesCurrRidge                  = EndPoints(currRidgeIm);
%locate the position
indexCurrEdges                  = find (edgesCurrRidge);
%create two matrices, each with one of the end points
borders411                       = zeros(size(currRidgeIm));
borders421                       = zeros(size(currRidgeIm));
borders411(indexCurrEdges(1))    = 1;
borders421(indexCurrEdges(2))    = 1;

[rB41,cB41]                     = find(borders41);
[rB42,cB42]                     = find(borders42);
[rB411,cB411]                     = find(borders411);
[rB421,cB421]                     = find(borders421);

dist11_22 = sqrt((rB41-rB411)^2+(cB41-cB411)^2)+sqrt((rB42-rB421)^2+(cB42-cB421)^2);
dist12_21 = sqrt((rB41-rB421)^2+(cB41-cB421)^2)+sqrt((rB42-rB411)^2+(cB42-cB411)^2);
%%
if dist12_21<dist11_22
    borders41=borders421;
    borders42=borders411;
else
    borders41=borders411;
    borders42=borders421;
    
end

%%




qqq1=bwdist(borders41);
qqq2=bwdist(borders42);


% figure(5)
% subplot(121)
% surfdat(qqq1+100*extendedEdge11)
% subplot(122)
% surfdat(qqq2+100*extendedEdge22)

%%

for counterEdge =1:numExEdges11
    www=qqq1.*(extendedEdge11==counterEdge);
    www2=(find(www));
    www3=sort(www(www2));
    minExtendedEdge11(counterEdge)=www3(1);
end
[a,b1]=min(minExtendedEdge11);
for counterEdge =1:numExEdges22
    www=qqq2.*(extendedEdge22==counterEdge);
    www2=(find(www));
    www3=sort(www(www2));
    minExtendedEdge22(counterEdge)=www3(1);

end
[a,b2]=min(minExtendedEdge22);
extendedEdgeCombined5           = (imdilate((extendedEdge11==b1)|(extendedEdge22==b2),ones(5))|imdilate(currRidgeIm,ones(4)));  %[1 1; 1 1]));

% surfdat(extendedEdgeCombined5)

twoRegions                      = bwlabel(1-bwmorph(extendedEdgeCombined5,'fill'),4);
% figure(3)
% surfdat(twoRegions+4*currRidgeIm)
% drawnow
%pause(0.05)
bordersEdge1                    = currRidgeIm.*imdilate(borders41,dilFactor3{5});
bordersEdge2                    = currRidgeIm.*imdilate(borders42,dilFactor3{5});

%%

% surfdat(qqq.*extendedEdge111)


%%
% extendedEdge111                 = unique(extendedEdge11.*imdilate(borders41,ones(10)));
% extendedEdge222                 = unique(extendedEdge22.*imdilate(borders42,ones(10)));
% extendedEdgeCombined1           = ismember(extendedEdge11,extendedEdge111(2:end));
% extendedEdgeCombined2           = ismember(extendedEdge22,extendedEdge222(2:end));
% extendedEdgeCombined5           = (imdilate(extendedEdgeCombined1|extendedEdgeCombined2,ones(3))|imdilate(currRidgeIm,[1 1; 1 1]));
% twoRegions                      = bwlabel(1-extendedEdgeCombined5,4);
% surfdat(twoRegions+4*currRidgeIm)


% bordersEdge1                    = currRidgeIm.*imdilate(borders41,dilFactor3{5});
% bordersEdge2                    = currRidgeIm.*imdilate(borders42,dilFactor3{5});


% %%
% extendedEdgeCombined1           = (extendedEdge22|extendedEdge11);
% extendedEdgeCombined2           = bwlabel(extendedEdgeCombined1.*(1-imdilate(currRidgeIm,[1 1;1 1])));
% extendedEdgeCombined31          = unique(extendedEdgeCombined2.*imdilate(borders41|borders42,ones(5)));
% extendedEdgeCombined3           = unique(extendedEdgeCombined2.*imdilate(borders41|borders42,ones(5)));
% extendedEdgeCombined4           = ismember(extendedEdgeCombined2,extendedEdgeCombined3(2:end));
% %extendedEdgeCombined2            = imdilate(extendedEdgeCombined1|currRidgeIm,[1 1;1 1]);
% %extendedEdgeCombined1            = imdilate(extendedEdge22|extendedEdge11|currRidgeIm,[1 1;1 1]);
% extendedEdgeCombined5           = (imdilate(extendedEdgeCombined4,[0 1 0;1 1 1;0 1 1])|imdilate(currRidgeIm,[1 1; 1 1]));
% 
% twoRegions                      = bwlabel(1-extendedEdgeCombined5,4);
% figure(1)
% surfdat(twoRegions+4*currRidgeIm)

%% Previous attempt to generate two regions by watershed of the distance transform
% %disp(sizeThicLongRidges(currRidge).Area)
%  qqq=-bwdist(currRidgeIm);
% minDistV= min(qqq(:))-1;
% %%
% if (abs(sizeThicLongRidges(currRidge).Orientation)<45)
%     qqq(1,:)                    = minDistV;
%     qqq(end,:)                  = minDistV;
% else
%     qqq(:,1)                    = minDistV;
%     qqq(:,end)                  = minDistV;    
% end

%%
% surfdat(imfilter(qqq,ones(13),'replicate'))

%%
%  twoRegions = (watershed(imfilter(qqq,ones(38))));
%  figure(3)
%  surfdat(twoRegions)
 
% surfdat(twoRegions+4*currRidgeIm)
% %%
% if (max(twoRegions(:))>2)    
% %%    
%     twoRegionsLines             = (twoRegions==0);
%     %test for holes in the watershed
%     twoRegionsLinesHoles        = imfill(twoRegionsLines,'holes');
%     if any (twoRegionsLines(:)~=twoRegionsLinesHoles(:))  
%         remHole                 = imdilate(-twoRegionsLines+twoRegionsLinesHoles,ones(3));
%         twoRegionsLines         = currRidgeIm|twoRegionsLines&(~remHole);
%         
%     end
%     twoRegionsLines             = bwmorph(twoRegionsLines,'spur',inf);
%     [currBranches,numBranches]  = bwlabel(twoRegionsLines-BranchPoints(twoRegionsLines));
%     qqq2=regionprops(currBranches,'Area');
%     %surfdat(currBranches);
% %%    
%     while (numBranches>1)
% %%        
%          %currBranches=ismember(currBranches,find([qqq2.Area]~=min([qqq2.Area])));
%          branchesWithMinSize    = find([qqq2.Area]==min([qqq2.Area]));
%          remBranch              = ismember(currBranches,branchesWithMinSize(1));
%          twoRegionsLines        = twoRegionsLines - remBranch;
%          [currBranches,numBranches]   = bwlabel(twoRegionsLines-BranchPoints(twoRegionsLines));
%          qqq2=regionprops(currBranches,'Area');
%          %numBranches = numBranches -1;
%          %surfdat(twoRegionsLines)
% %%
%     end
% %%
%     twoRegions = bwlabel(1-bwmorph(currBranches,'diag'));%(watershed(currBranches));
% %%    
% end

%%



% for the stop line, add 90 degrees to get a perpendicular line and repeat without spur


edgesCurrRidge                  = EndPoints(currRidgeIm);
%locate the position
indexCurrEdges                  = find (edgesCurrRidge);
%create two matrices, each with one of the end points
borders41                       = zeros(size(currRidgeIm));
borders42                       = zeros(size(currRidgeIm));
borders41(indexCurrEdges(1))    = 1;
borders42(indexCurrEdges(2))    = 1;
%
%trim the ridge to the last 5-10 points to get the orientation of the ridge CLOSE to the endpoint
bordersEdge1                    = currRidgeIm.*imdilate(borders41,dilFactor3{5});
bordersEdge2                    = currRidgeIm.*imdilate(borders42,dilFactor3{5});
angle1                          = regionprops(bwlabel(bordersEdge1),'orientation');
angle2                          = regionprops(bwlabel(bordersEdge2),'orientation');





EndPoint1_perpend               = strel('line', 80, 90+angle1(1).Orientation );
EndPoint2_perpend               = strel('line', 80, 90+angle2(1).Orientation);



[borders2]          =  edge(currData2,'canny',[],2);
[borders2L]         = bwlabel(borders2);
%overlapBord             = (borders2L.*(imdilate(currRidgeIm,ones(5))));
%keepBorders         = setdiff((1:numBordL),unique(overlapBord));
%borders22           = ismember(borders2L,(1:numBordL));

borders22               = (borders2L.*(1-imdilate(currRidgeIm,dilFactor2)));
borders22L          = bwlabel(bwmorph(borders22,'spur',1));
tempAreaBor         = regionprops(borders22L,'Area');
borders3            = ismember(borders22L,find([tempAreaBor.Area]>8));
borders3            = borders3|borders1;

%%
% figure(23)
% subplot(121)
% surfdat(2*currRidgeIm+((twoRegions==1).*(borders3+fRidges2D).*(1-currRidgeIm))>0)
% subplot(122)
% surfdat(2*currRidgeIm+((twoRegions==2).*(borders3+fRidges2D).*(1-currRidgeIm))>0)


%%
distancesSide1      = bwdist((twoRegions==1).*(borders3+fRidges2D).*(1-currRidgeIm));
distancesSide2      = bwdist((twoRegions==2).*(borders3+fRidges2D).*(1-currRidgeIm));
%%
% figure(2)
% subplot(121)
% surfdat(distancesSide1)
% subplot(122)
% surfdat(distancesSide2)
%%
side1Dist           = currRidgeIm.*distancesSide1;
side2Dist           = currRidgeIm.*distancesSide2;

%%
% dilFactor3{50}=0;
% for k=50:-1:1;
%     dilFactor3{k}   =strel('disk',k,0);
%     %dilFactorSize(k)=1+k*2;%size(getnhood(dilFactor3{k}),1);
% end
%%
borders4            = zeros(size(borders3));
bordersSide1        = zeros(size(borders3));
bordersSide2        = zeros(size(borders3));
indexCurrRidge      = find (side2Dist);
[indexCurrRidgeR,indexCurrRidgeC]      = find (side2Dist);
side1DistIn         = round(side1Dist(indexCurrRidge));
side2DistIn         = round(side2Dist(indexCurrRidge));

fBoundaryVal        = mean(side1DistIn+side2DistIn);

%%
% for k=2:6:sizeThicLongRidges(currRidge).Area
%     disp(k)
%     borders4(indexCurrRidge(k))     = 1;
%     bordersSide1 = bordersSide1 | imdilate(borders4,dilFactor3{round(side1DistIn((k)))});
%     bordersSide2 = bordersSide2 | imdilate(borders4,dilFactor3{round(side2DistIn((k)))});
% 
%   %  bordersSide1(indexCurrRidgeR(k),indexCurrRidgeC(k)) = bordersSide1(indexCurrRidgeR(k),indexCurrRidgeC(k)) | getnhood(dilFactor3{side1DistIn((k))});
%   %  bordersSide2(indexCurrRidgeR(k),indexCurrRidgeC(k)) = bordersSide2(indexCurrRidgeR(k),indexCurrRidgeC(k)) | getnhood(dilFactor3{side2DistIn((k))});
%     
%    
%     borders4(indexCurrRidge(k))     = 0;
% end
%%


[rowsCurr,colsCurr ]= size(currRidgeIm);
for k=2:2:sizeThicLongRidges(currRidge).Area
    %disp(k)
    %borders4(indexCurrRidge(k))     = 1;
    %  bordersSide1 = bordersSide1 | imdilate(borders4,dilFactor3{round(side1DistIn((k)))});
    %  bordersSide2 = bordersSide2 | imdilate(borders4,dilFactor3{round(side2DistIn((k)))});
    try
        kDil = floor(size(getnhood(dilFactor3{min(side1DistIn(k),60)}),1)/2);
    catch

        qqq=1; %#ok<NASGU>
    end
    rDimOver = indexCurrRidgeR(k)-kDil:indexCurrRidgeR(k)+kDil;
    cDimOver = indexCurrRidgeC(k)-kDil:indexCurrRidgeC(k)+kDil;

    kernelToOverlap = getnhood(dilFactor3{min(side1DistIn(k),60)});
    kernelToOverlap = kernelToOverlap (find(ismember(rDimOver,1:rowsCurr)),find(ismember(cDimOver,1:colsCurr)));


    rDimOver = rDimOver (find(ismember(rDimOver,1:rowsCurr)));
    cDimOver = cDimOver (find(ismember(cDimOver,1:colsCurr)));


    try
        bordersSide1(rDimOver,cDimOver) = bordersSide1(rDimOver,cDimOver) |kernelToOverlap ;

    catch
        qqq=1; %#ok<NASGU>
    end
    try
        kDil = floor(size(getnhood(dilFactor3{min(side2DistIn(k),60)}),1)/2);
    catch
        qqq=1; %#ok<NASGU>
    end
    rDimOver = indexCurrRidgeR(k)-kDil:indexCurrRidgeR(k)+kDil;
    cDimOver = indexCurrRidgeC(k)-kDil:indexCurrRidgeC(k)+kDil;

    try
        kernelToOverlap = getnhood(dilFactor3{min(side2DistIn(k),60)});
    catch
        qqq=1; %#ok<NASGU>
    end
    kernelToOverlap = kernelToOverlap (find(ismember(rDimOver,1:rowsCurr)),find(ismember(cDimOver,1:colsCurr)));


    rDimOver = rDimOver (find(ismember(rDimOver,1:rowsCurr)));
    cDimOver = cDimOver (find(ismember(cDimOver,1:colsCurr)));

    try
        bordersSide2(rDimOver,cDimOver) = bordersSide2(rDimOver,cDimOver) |kernelToOverlap ;

    catch
        qqq=1; %#ok<NASGU>
    end
    borders4(indexCurrRidge(k))     = 0;
end



%% calculate the borders according to the dilated endpoint
bordersEdge11                   = imdilate(borders41,EndPoint1_perpend);
bordersEdge22                   = imdilate(borders42,EndPoint2_perpend);
bordersEdgeCombined             = imdilate(bordersEdge22|bordersEdge11,[1 1;1 1]);


borders2SidesFull               = bwlabel((bordersSide1.*(twoRegions==1) + bordersSide2.*(twoRegions==2)+currRidgeImDil ).*(1-bordersEdgeCombined));
sizeFullBorders                 = regionprops(borders2SidesFull,'Area');

[tt1,tt2]                       = sort([sizeFullBorders.Area]);

borders2Sides                   = ismember(borders2SidesFull,tt2(end))+(imdilate(currRidgeIm,ones(3)));

finalBoundary                   = imclose(borders2Sides>0,ones(3));
%%
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [finalBoundary]   = findVessBoundary3(dilFactor3,fRidges2DthicLong_L,ridgeWidth)
%This function is for those ridges that are thin and the estimation of thickness is uniform and obtained
%from the depth of the scale space estimation

%%
% dilFactor3{50}=0;
% for k=50:-1:1;
%     dilFactor3{k}   =strel('disk',k,0);
%     %dilFactorSize(k)=1+k*2;%size(getnhood(dilFactor3{k}),1);
% end

side1DistIn         = round(ridgeWidth);

currRidgeIm             = bwmorph(fRidges2DthicLong_L,'spur',1+ceil(side1DistIn/2));


bordersSide1        = zeros(size(currRidgeIm));

%%
sizeCurrRidge       = sum(currRidgeIm(:));

[indexCurrRidgeR,indexCurrRidgeC]      = find (currRidgeIm);

[rowsCurr,colsCurr ]= size(currRidgeIm);
for k=2:1:sizeCurrRidge
    try
        kDil = floor(size(getnhood(dilFactor3{side1DistIn}),1)/2);
    catch
        qqq=1; %#ok<NASGU>
    end
    %find dimensions of the region to be overlapped by the dilated disk over the ridge
    rDimOver = indexCurrRidgeR(k)-kDil:indexCurrRidgeR(k)+kDil;
    cDimOver = indexCurrRidgeC(k)-kDil:indexCurrRidgeC(k)+kDil;
    
    % the kernel has been previously defined, but limit dimensions to [1 rows]x[1 cols]
    kernelToOverlap = getnhood(dilFactor3{side1DistIn});
    kernelToOverlap = kernelToOverlap (find(ismember(rDimOver,1:rowsCurr)),find(ismember(cDimOver,1:colsCurr)));

    % new dimensions
    rDimOver = rDimOver (find(ismember(rDimOver,1:rowsCurr)));
    cDimOver = cDimOver (find(ismember(cDimOver,1:colsCurr)));

    % logical or between the existing region and the kernel
    try
        bordersSide1(rDimOver,cDimOver) = bordersSide1(rDimOver,cDimOver) |kernelToOverlap ;
        
%          figure(5)
%          surfdat(bordersSide1+fRidges2DthicLong_L)
%          drawnow
%          pause(0.1)
    catch
        qqq=1; %#ok<NASGU>
    end
end
finalBoundary                   = bordersSide1;

%% Limit side semicircular edges by cutting with a perpendicular line to the end point

% clear angle1 angle2;
% 
% edgesCurrRidge                  = EndPoints(currRidgeIm);
% indexCurrEdges                  = find (edgesCurrRidge);
% 
% borders41                       = zeros(size(currRidgeIm));
% try
% borders41(indexCurrEdges(1))    = 1;
% catch
%     qqq=1;
% end
% bordersEdge1                    = currRidgeIm.*imdilate(borders41,dilFactor3{5});
% angle1                          = regionprops(bwlabel(bordersEdge1),'orientation');
% SE1                             = strel('line', 80,90+angle1(1).Orientation );
% bordersEdge11                   = imdilate(borders41,SE1);
% 
% borders42                       = zeros(size(currRidgeIm));
% borders42(indexCurrEdges(2))    = 1;
% 
% bordersEdge2                    = currRidgeIm.*imdilate(borders42,dilFactor3{5});
% angle2                          = regionprops(bwlabel(bordersEdge2),'orientation');
% SE2                             = strel('line', 80, 90+angle2(1).Orientation);
% bordersEdge22                   = imdilate(borders42,SE2);
% 
% 
% 
% bordersEdgeCombined             = imdilate(bordersEdge22|bordersEdge11,[1 1;1 1]);
% 
% borders2SidesFull               = bwlabel((bordersSide1 ).*(1-bordersEdgeCombined));
% sizeFullBorders                 = regionprops(borders2SidesFull,'Area');
% 
% [tt1,tt2]                       = sort([sizeFullBorders.Area]);
% 
% borders2Sides                   = ismember(borders2SidesFull,tt2(end-1:end))+(imdilate(currRidgeIm,ones(3)));
% 
% finalBoundary                   = imclose(borders2Sides>0,ones(3));
%%
end
