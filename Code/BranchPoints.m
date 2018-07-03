function [BranchPoints1,BranchPoints2,numBranchPoints]=BranchPoints(dataIn)

BranchKernel1                   = [0 1 0; 1 1 0; 0 0 1];    %Corner
BranchKernel2                   = [0 1 0; 0 1 0; 1 0 1];    %Y
BranchKernel3                   = [1 0 0; 0 1 0; 1 0 1];    %diagonal
%BranchKernel4                  = [1 1 1; 0 1 0; 0 1 0];    %T

dataIn=(padData(dataIn,1,[],0));

BranchPoints1=zeros(size(dataIn));

for k=0:3
    BranchPoints1               = (BranchPoints1|bwhitmiss(dataIn,rot90(BranchKernel1,k)));
    BranchPoints1               = (BranchPoints1|bwhitmiss(dataIn,rot90(BranchKernel2,k)));
    BranchPoints1               = (BranchPoints1|bwhitmiss(dataIn,rot90(BranchKernel3,k)));
end

BranchPoints1                   = BranchPoints1(2:end-1,2:end-1);



if nargout>=2
    BranchPoints2               = BranchPoints1;
    BranchPoints2(1:end-1,:)    = (BranchPoints1(1:end-1,:)|BranchPoints1(2:end,:));
    BranchPoints2(2:end,:)      = (BranchPoints1(1:end-1,:)|BranchPoints2(2:end,:));
    BranchPoints2(:,1:end-1)    = (BranchPoints2(:,1:end-1)|BranchPoints1(:,2:end));
    BranchPoints2(:,2:end)      = (BranchPoints1(:,1:end-1)|BranchPoints2(:,2:end));
    
    numBranchPoints             = sum(find(BranchPoints1)>0);
    
end

%BranchPoints=(BranchPoints1|BranchPoints2|BranchPoints3);