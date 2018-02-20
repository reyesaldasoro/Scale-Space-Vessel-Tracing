function [EndPoints1,EndPoints2,numEndPoints]=EndPoints(dataIn)
% [EndPoints1,EndPoints2,numEndPoints]=EndPoints(dataIn)
%-------------------------------------------------------------------------------------
%------  Author :   Constantino Carlos Reyes-Aldasoro                       ----------
%------             The University of Sheffield                             ----------
%------             http://carlos-reyes.staff.shef.ac.uk                    ----------
%------  22 Oct 2009                                       ---------------------------
%-------------------------------------------------------------------------------------


EndKernel1                   = [ 0 0 0; -1 1 -1; -1 -1 -1];    %End points
EndKernelEx                  = [ 1 0 1; -1 1 -1; -1 -1 -1];    %End points

%EndKernel2                   = [0 1 0; 0 1 0; 1 0 1];    %Y
%EndKernel3                   = [1 0 0; 0 1 0; 1 0 1];    %diagonal
%EndKernel4                  = [1 1 1; 0 1 0; 0 1 0];    %T

dataIn=(padData(dataIn,1,[],0));

EndPoints1                      =zeros(size(dataIn));
EndPointsEx                     = zeros(size(dataIn));
%%
for k=0:3
    EndPoints1               = (EndPoints1|bwhitmiss(dataIn,rot90(EndKernel1,k)));
   % EndPoints1               = (EndPoints1|bwhitmiss(dataIn,rot90(EndKernel2,k)));
   % EndPoints1               = (EndPoints1|bwhitmiss(dataIn,rot90(EndKernel3,k)));
end

EndPoints1                   = EndPoints1(2:end-1,2:end-1);
%%


for k=0:3
    EndPointsEx               = (EndPointsEx|bwhitmiss(dataIn,rot90(EndKernelEx,k)));
   % EndPoints1               = (EndPoints1|bwhitmiss(dataIn,rot90(EndKernel2,k)));
   % EndPoints1               = (EndPoints1|bwhitmiss(dataIn,rot90(EndKernel3,k)));
end

EndPointsEx                   = EndPointsEx(2:end-1,2:end-1);

%%
EndPoints1 = EndPoints1 &(~EndPointsEx);

%%
if nargout>=2
    EndPoints2               = EndPoints1;
    EndPoints2(1:end-1,:)    = (EndPoints1(1:end-1,:)|EndPoints1(2:end,:));
    EndPoints2(2:end,:)      = (EndPoints1(1:end-1,:)|EndPoints2(2:end,:));
    EndPoints2(:,1:end-1)    = (EndPoints2(:,1:end-1)|EndPoints1(:,2:end));
    EndPoints2(:,2:end)      = (EndPoints1(:,1:end-1)|EndPoints2(:,2:end));
    
    numEndPoints             = sum(find(EndPoints1)>0);
    
end

%EndPoints=(EndPoints1|EndPoints2|EndPoints3);