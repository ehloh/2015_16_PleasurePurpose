function [ V ] = b01(x, modelinput)
%
%       x:                    Free parameters       (row=parampoint, col=parameter)
%       modelinput:     {[]    data   fPar   col}
%        V:                   Value associated with Opt 1 & Opt 2
%                                   (row parampoint, col=trialnum, 3rd dimension=Option)
%
% ---------------------------------------------------------------------------------------------------------------------------------

data=modelinput{2}; fPar=modelinput{3}; col=modelinput{4}; nPP=size(x,1); nTrials=size(data,1);  % No. of parameter points (unique combinations of parameter values for all free parameters covered by grid)
pl1= shiftdim(data(:, col.PL(1)),1);  pp1=shiftdim(data(:, col.PP(1)),1) ; pl2= shiftdim(data(:, col.PL(2)),1); pp2=   shiftdim(data(:, col.PP(2)),1);

%%  Value calculation 
 
% Value just considers PL scores  
V(:,:,1) = repmat(pl1, nPP,1);
V(:,:,2) = repmat(pl2, nPP,1);

end


