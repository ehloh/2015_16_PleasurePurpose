function [ V ] = b01(x, modelinput)
%
%       x:                    Free parameters       (row=parampoint, col=parameter)
%       modelinput:     {[]    data   fPar   col}
%        V:                   Value associated with Opt 1 & Opt 2
%                                   (row parampoint, col=trialnum, 3rd dimension=Option)
%
% ---------------------------------------------------------------------------------------------------------------------------------

aw= 4/(1+exp(-x(:,2)))-2;

%
data=modelinput{2}; fPar=modelinput{3}; col=modelinput{4}; nPP=size(x,1); nTrials=size(data,1);  % No. of parameter points (unique combinations of parameter values for all free parameters covered by grid)
pl1= shiftdim(data(:, col.PL(1)),1);  pp1=shiftdim(data(:, col.PP(1)),1) ; pl2= shiftdim(data(:, col.PL(2)),1); pp2=   shiftdim(data(:, col.PP(2)),1);

%%  Value calculation 

plcf = shiftdim(data(:, col.PLcf,1),1);  
ppcf = shiftdim(data(:, col.PPcf,1),1);  


keyboard

% Determine weight of pleasure score 
%           Note: scale of weight should be <1 (for softmax, given constraints of aw)
plweight = repmat(aw, 1, nTrials).*repmat(0.1*plcf./ppcf, nPP,1);  %  

% Scale pleasure scores 
pl1= repmat(pl1, nPP,1).*plweight; 
pl2= repmat(pl2, nPP,1).*plweight; 



% % Linearly re-scale PL and PP scores to maintain (rough) 1-10 range 
% scoremax= max([pl1(:);pl2(:); pp1(:); pp2(:)]); 
% pl1= pl1.*(10./ scoremax); 
% pl2= pl2.*(10./ scoremax); 
% pp1= pp1.*(10./ scoremax); 
% pp2= pp2.*(10./ scoremax); 

%%

keyboard

% Linearly combine altered attribute scores 
V(:,:,1) = pl1+pp1; 
V(:,:,2) = pl2+pp2; 

end


