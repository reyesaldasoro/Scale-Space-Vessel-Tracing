function dataOut=removeImageOutliers(dataIn,sizeNeig)

%some images have salt-and-noise noise on the edges (extreme values not related to the data).
%Remove them by replacing with the median all those values that exceed the 0.2% values of the mean
% (0-0.2) and (0.998-1)

% usual dimension check
[rows,cols,levs]                    = size(dataIn);

if ~exist('sizeNeig','var')
    sizeNeig =3;
end

if levs>1
    for counterLevs=1:levs
        currLevel                           = dataIn(:,:,counterLevs);
        %currLevel0                          = imfilter(currLevel,gaussF(5,5,1));
        avcurrLevel                         = repmat(mean2(currLevel),[rows cols]);
        q3                                  = numel(currLevel);
        [q1]                                = (sort(currLevel(:)));
        min002                              = q1(round(q3*0.001));
        max998                              = q1(round(q3*0.999));
        %lowCurrLevel                        = repmat(q1(round(q3*0.01)),[rows cols]);
        %highCurrLevel                       = repmat(q1(round(q3*0.99)),[rows cols]);
        currLevel2                          = medfilt2(padData(currLevel,1),[sizeNeig sizeNeig]);
        currLevel2                          = currLevel2(2:end-1,2:end-1);
        currLevel(currLevel<=min002)        = max(avcurrLevel(currLevel<=min002),currLevel2(currLevel<=min002));
        currLevel(currLevel>=max998)        = min(avcurrLevel(currLevel>=max998),currLevel2(currLevel>=max998));%currLevel2(currLevel>=max998);

        %currLevel(currLevel<=min002)        =
        %currLevel(currLevel<=min002)        =
        dataOut(:,:,counterLevs)            = currLevel; %#ok<AGROW>
    end
else

    dataIn0                             = dataIn(3:end-2,3:end-2);
    q3                                  = numel(dataIn0);
    [q1]                                = (sort(dataIn0(:)));
    min002                              = q1(round(q3*0.005));
    max998                              = q1(round(q3*0.995));
    %%
    %Create a new copy of the data, with two extra rows and columns as the median filter pads with zeros
    %copy the SECOND line to the first to avoid the noise duplicating

    dataIn3                             = zeros(rows+2,cols+2);
    dataIn3(2:end-1,2:end-1)            = dataIn;

    dataIn3(1,2:end-1)                  = dataIn(3,:);
    dataIn3(end,2:end-1)                = dataIn(end-2,:);
    dataIn3(2:end-1,1)                  = dataIn(:,3);
    dataIn3(2:end-1,end)                = dataIn(:,end-2);
    dataIn3(1,1)                        = dataIn(3,3);
    dataIn3(1,end)                      = dataIn(3,end);
    dataIn3(end,1)                      = dataIn(end,3);
    dataIn3(end,end)                    = dataIn(end-1,end-1);


    %proceed to filter
    dataIn3                             = medfilt2(dataIn3,[sizeNeig sizeNeig]);

    %remove edges
    dataIn4                             = dataIn3(2:end-1,2:end-1);
    %% Remove outliers
    dataOut                             = dataIn;
    dataOut(dataIn<min002)              = dataIn4(dataIn<min002);
    dataOut(dataIn>max998)              = dataIn4(dataIn>max998);

    %% plus, the edges will be replaced by medians

    dataOut(1,:)                        = dataIn4(1,:);
    dataOut(end,:)                      = dataIn4(end,:);
    dataOut(:,1)                        = dataIn4(:,1);
    dataOut(:,end)                      = dataIn4(:,end);


end

