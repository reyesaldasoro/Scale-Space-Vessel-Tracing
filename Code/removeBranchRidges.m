function [fRidges,numRidges] = removeBranchRidges(fRidges)


for counterPixels=1:2:11
    %figure(1); imagesc(sum(fRidges,3))
    fRidgesOut                                      = fRidges ;
    ridges2D                                        = sum(fRidgesOut,3);

    %detect branch points
    [BranchPoints1,BranchPoints2,numBranchPoints]   = BranchPoints(ridges2D);
    [rr2,cc2]                                       = find(BranchPoints1);

    %remove branch points from 3D structure
    for counterBP =1:numBranchPoints
        fRidgesOut(rr2(counterBP),cc2(counterBP),:)    = 0; %#ok<AGROW>
    end

    %label again and remove objects that are smaller than X pixels

    fRidgesL                                        = bwlabeln(fRidgesOut>0);
    propsRidges                                     = regionprops(fRidgesL,'Area');
    fRidgesToR                                      = ismember(fRidgesL,find([propsRidges.Area] <= counterPixels));
    fRidges                                         = bwlabeln((fRidges>0) - (fRidgesToR>0));
    %figure(2); imagesc(sum(fRidges,3))
    %beep;
end

ridges2D                                        = sum(fRidges,3);

%detect branch points
[BranchPoints1,BranchPoints2,numBranchPoints]   = BranchPoints(ridges2D);
[rr2,cc2]                                       = find(BranchPoints1);

%remove branch points from 3D structure
for counterBP =1:numBranchPoints
    fRidges(rr2(counterBP),cc2(counterBP),:)    = 0; %#ok<AGROW>
end
fRidgesL                                        = bwlabeln(fRidges>0);
propsRidges                                     = regionprops(fRidgesL,'Area');
fRidgesToKeep                                   = ismember(fRidgesL,find([propsRidges.Area] >3));
[fRidges,numRidges]                             = bwlabeln((fRidgesToKeep>0) );

