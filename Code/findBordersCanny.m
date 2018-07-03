[rows,cols,levs] =  size(dataIn);

numRidges = size(fStats,1);



fRidges2Dthin=zeros(rows-2,cols-2);
fRidges2Dthic=zeros(rows-2,cols-2);

for k=1:numRidges
    tempRidge       = (fRidges==k);
    %tempRidgeV      = posInRidge.*(tempRidge);
    tempRidge2D     = sum(tempRidge,3);
    %tempLevel = mean(ttt(find(ttt)));
    if fStats(k,3)>6
        fRidges2Dthic =fRidges2Dthic | tempRidge2D;
    else
        fRidges2Dthin =fRidges2Dthin | tempRidge2D;
    end
    
    
end

fRidges2Dthin           = padData(fRidges2Dthin,1);
fRidges2Dthic           = padData(fRidges2Dthic,1);

fRidges2D               = (fRidges2Dthin|fRidges2Dthic);

fRidges2Dthic_L         = bwlabel(fRidges2Dthic);
fRidges2Dthin_L         = bwlabel(fRidges2Dthin);

sizeThinRidges          = regionprops(fRidges2Dthin_L,'Area','Centroid','BoundingBox','orientation');
sizeThicRidges          = regionprops(fRidges2Dthic_L,'Area','Centroid','BoundingBox','orientation');

%%
fRidges2DthicLong_L       = bwlabel(ismember(fRidges2Dthic_L,find([sizeThicRidges.Area]>40)));
sizeThicLongRidges      = regionprops(fRidges2DthicLong_L,'Area','Centroid','BoundingBox','orientation');

numThicR                = numel(sizeThicRidges);
numThinR                = numel(sizeThinRidges);
numThicLongR            = numel(sizeThicLongRidges);

%%

dilFactor0  = 45;
dilFactor   =strel('disk',dilFactor0,0);
dilFactor2  =strel('disk',5,0);


%%

currRidge =30;

disp(sizeThicLongRidges(currRidge).Orientation)
imagesc (fRidges2DthicLong_L==currRidge);



%%
currRidge               = 30;

currRowsD = max(1,ceil(sizeThicLongRidges(currRidge).BoundingBox(2))-dilFactor0):min(rows,floor(sizeThicLongRidges(currRidge).BoundingBox(2))+(sizeThicLongRidges(currRidge).BoundingBox(4))+dilFactor0);
currColsD = max(1,ceil(sizeThicLongRidges(currRidge).BoundingBox(1))-dilFactor0):min(cols,floor(sizeThicLongRidges(currRidge).BoundingBox(1))+(sizeThicLongRidges(currRidge).BoundingBox(3))+dilFactor0);

currRows = (ceil(sizeThicLongRidges(currRidge).BoundingBox(2))):(floor(sizeThicLongRidges(currRidge).BoundingBox(2))+(sizeThicLongRidges(currRidge).BoundingBox(4)));
currCols = (ceil(sizeThicLongRidges(currRidge).BoundingBox(1))):(floor(sizeThicLongRidges(currRidge).BoundingBox(1))+(sizeThicLongRidges(currRidge).BoundingBox(3)));


%disp(sizeThicLongRidges(currRidge).Area)
currRidgeIm             = (fRidges2DthicLong_L(currRowsD,currColsD)==currRidge);
qqq=-bwdist(currRidgeIm);

minDistV= min(qqq(:))-1;
if (abs(sizeThicLongRidges(currRidge).Orientation)<45)
    qqq(1,:)                    = minDistV;
    qqq(end,:)                  = minDistV;
    %qqq(1:floor(end/3)-1,1)     = minDistV;
    %qqq(1:floor(end/3)-1,end)   = minDistV;
    %qqq(ceil(end/3)+1:end,1)    = minDistV;
    %qqq(ceil(end/3)+1:end,end)  = minDistV;
else
    qqq(:,1)                    = minDistV;
    qqq(:,end)                  = minDistV;    
    %qqq(1,1:floor(end/3)-1)     = minDistV;
    %qqq(end,1:floor(end/3)-1)   = minDistV;
    %qqq(1,ceil(end/3)+1:end)    = minDistV;
    %qqq(end,ceil(end/3)+1:end)  = minDistV;
end
twoRegions = (watershed(qqq));

if (max(twoRegions(:)>2))
    
    [currBranches,numBranches]   = bwlabel((twoRegions==0)-BranchPoints(twoRegions==0));
    qqq2=regionprops(currBranches,'Area');
    while (numBranches>1)
         currBranches=ismember(currBranches,find([qqq2.Area]~=min([qqq2.Area])));
         numBranches = numBranches -1;
    end
    twoRegions = (watershed(currBranches));
    
end

%%
currData                = dataInSm(currRowsD,currColsD).*imdilate(currRidgeIm,dilFactor);
currData2               = currData;
currData2(currData==0)  = 999;
minData2                = min(currData2(:))-1;
currData2(currData2==999)= minData2;

[borders2,thresC]       =  edge(currData2,'canny',[],2);
[borders2L,numBordL]    = bwlabel(borders2);
%overlapBord             = (borders2L.*(imdilate(currRidgeIm,ones(5))));
%keepBorders         = setdiff((1:numBordL),unique(overlapBord));
%borders22           = ismember(borders2L,(1:numBordL));

borders22               = (borders2L.*(1-imdilate(currRidgeIm,dilFactor2)));
borders22L          = bwlabel(bwmorph(borders22,'spur',1));
tempAreaBor         = regionprops(borders22L,'Area');
borders3            = ismember(borders22L,find([tempAreaBor.Area]>8));
%%
%imagesc((1-(borders3+fRidges2D(currRowsD,currColsD))).*dataInSm2(currRowsD,currColsD)+50*currRidgeIm)
%imagesc((twoRegions==1).*(borders3+fRidges2D(currRowsD,currColsD)).*(1-currRidgeIm))
distancesSide1      = bwdist((twoRegions==1).*(borders3+fRidges2D(currRowsD,currColsD)).*(1-currRidgeIm));
distancesSide2      = bwdist((twoRegions==2).*(borders3+fRidges2D(currRowsD,currColsD)).*(1-currRidgeIm));

%distancesEdge1      = bwdist(edgesCurrRidge));
%distancesEdge2      = bwdist((borders3+fRidges2D(currRowsD,currColsD)).*(edgesCurrRidge));
side1Dist           = currRidgeIm.*distancesSide1;
side2Dist           = currRidgeIm.*distancesSide2;





%imagesc(borders2+2*(fRidges2Dthic_L==currRidge))
%imagesc((1-(fRidges2Dthic_L==currRidge)).*currData)
%caxis([92 216])
%%
imagesc(side1Dist)
%%
for k=30:-1:1;
    dilFactor3{k} =strel('disk',k,0);
end
%%
borders4 = zeros(size(borders3));
bordersSide1 = zeros(size(borders3));
bordersSide2 = zeros(size(borders3));
indexCurrRidge = find (side2Dist);
side1DistIn         = side1Dist(indexCurrRidge);
side2DistIn         = side2Dist(indexCurrRidge);
%%
for k=1:sizeThicLongRidges(currRidge).Area
    disp(k)
    borders4(indexCurrRidge(k))     = 1;
    bordersSide1 = bordersSide1 | imdilate(borders4,dilFactor3{round(side1Dist(indexCurrRidge(k)))});
    bordersSide2 = bordersSide2 | imdilate(borders4,dilFactor3{round(side2Dist(indexCurrRidge(k)))});
    borders4(indexCurrRidge(k))     = 0;
    
    
end
%%


%%
edgesCurrRidge      = EndPoints(currRidgeIm);

indexCurrEdges      = find (edgesCurrRidge);
borders4 = zeros(size(borders3));
borders4(indexCurrEdges(1))     = 1;

bordersEdge1        = currRidgeIm.*imdilate(borders4,dilFactor3{5});
angle1              = regionprops(bwlabel(bordersEdge1),'orientation');
SE1                 = strel('line', 80,90+angle1.Orientation );
bordersEdge11       = imdilate(borders4,SE1);
borders4(indexCurrRidge(1))     = 0;
    
borders4(indexCurrEdges(2))     = 1;

bordersEdge2        = currRidgeIm.*imdilate(borders4,dilFactor3{5});
angle2              = regionprops(bwlabel(bordersEdge2),'orientation');
SE2                 = strel('line', 80, 90+angle2.Orientation);
bordersEdge22       = imdilate(borders4,SE2);

borders4(indexCurrRidge(2))     = 0;

bordersEdgeCombined = imdilate(bordersEdge22|bordersEdge11,[1 1;1 1]);
% figure(1)
%  surfdat(bordersEdge11)
% 
% figure(2)
%  surfdat(bordersEdge22)

borders2SidesFull = bwlabel((bordersSide1.*(twoRegions==1) + bordersSide2.*(twoRegions==2)).*(1-bordersEdgeCombined));
sizeFullBorders = regionprops(borders2SidesFull,'Area');

[tt1,tt2]=sort([sizeFullBorders.Area]);

borders2Sides = ismember(borders2SidesFull,tt2(end-1:end));

finalBoundary = imfill(borders2Sides+(imdilate(edgesCurrRidge,ones(6))).*(bordersEdgeCombined),'holes');

%%
edgesSide1          = side1Dist.*edgesCurrRidge;
edgesSide2          = side2Dist.*edgesCurrRidge;

edge1Side1          = imdilate(edgesCurrRidge,dilFactor3{round(side1Dist(indexCurrRidge(k)))});


%%



%%

surfdat(borders3+borders2Sides+fRidges2D(currRowsD,currColsD)+currRidgeIm*2)


%%
qq=zeros(64);
qq(32,32)=1;
%%
figure(3)
imagesc(imdilate(qq,strel('disk',25,0))-qq)
grid on
%%
sizeF               = 9;
filtg               = gaussF(sizeF,sizeF,1);
dataInSm            = double(imfilter(dataIn,filtg,'replicate'));

%fRidges2D           = padData(sum(fRidges>0,3),1);
%fRidges2D_L         = bwlabel(fRidges2D);
dataInSm2           = dataInSm.*(1-fRidges2Dthic) + ((15+ dataInSm).*fRidges2Dthic);
%%
borders1            = (watershed(bwdist(fRidges2Dthic))==0);
%[Lx,Ly]             = gradient(255-dataInSm);

%borders             = abs(Ly)+abs(Lx);
%surfdat(borders)
%%

[borders2,thresC]   =  edge(dataInSm(60:190,450:530),'canny',[0.096 0.24],2.48832);

figure(20)
imagesc(borders2+2*fRidges2Dthic(60:190,450:530))
%%
borders22           = (borders2.*(1-imdilate(fRidges2Dthic,ones(3))));
borders22L          = bwlabel(bwmorph(borders22,'spur',1));
tempAreaBor         = regionprops(borders22L,'Area');
borders3            = ismember(borders22L,find([tempAreaBor.Area]>6));

figure(3);
imagesc(2*borders3+3*fRidges2Dthic)
figure(2);
imagesc(2*borders1+3*fRidges2Dthic)

%%

figure(1)
imagesc(3*fRidges2Dthic+(1*borders1)+(2*borders3))
%%
figure(5)
imagesc(dataInSm2)
%L=watershed(255-borders);
%%
%clf
surfdat(dataInSm.*(1-fRidges2Dthic))

%%

numRidges = size(fStats,1);

posInRidge=repmat(permute(linspace(1,10,10),[1 3 2]),[rows-2 cols-2 1]);
%%

fRidges2Dthin=zeros(rows-2,cols-2);
fRidges2Dthic=zeros(rows-2,cols-2);

for k=1:numRidges
    tempRidge       = (fRidges==k);
    %tempRidgeV      = posInRidge.*(tempRidge);
    tempRidge2D     = sum(tempRidge,3);
    tempLevel = mean(ttt(find(ttt)));
    if fStats(k,3)>7
        fRidges2Dthic =fRidges2Dthic | tempRidge2D;
    else
        fRidges2Dthin =fRidges2Dthin | tempRidge2D;
    end
    
    
end

fRidges2Dthin=padData(fRidges2Dthin,1);
fRidges2Dthic=padData(fRidges2Dthic,1);



%%
figure(4)
clf
surfdat(fRidges2Dthic_L==5)


%%
figure(3)
imagesc(dataOut2)


    