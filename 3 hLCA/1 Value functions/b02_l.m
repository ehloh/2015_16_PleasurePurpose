function [ V ] = b01(x, modelinput)
%
%       x:                    Free parameters       (row=parampoint, col=parameter)
%       modelinput:     {[]    data   fPar   col}
%        V:                   Value associated with Opt 1 & Opt 2
%                                   (row parampoint, col=triplwnum, 3rd dimension=Option)
%
% ---------------------------------------------------------------------------------------------------------------------------------

plw=1./(2+2.*exp(-x(:,2)));


data=modelinput{2}; fPar=modelinput{3}; col=modelinput{4}; nPP=size(x,1);   nTrials=size(data,1); % No. of parameter points (unique combinations of parameter values for all free parameters covered by grid)
pl1= shiftdim(data(:, col.PL(1)),1);  pp1=shiftdim(data(:, col.PP(1)),1) ; pl2= shiftdim(data(:, col.PL(2)),1); pp2=   shiftdim(data(:, col.PP(2)),1);

%%  Value calculation 

V(:,:,1) = repmat(pl1, nPP,1).*repmat(plw, 1,nTrials)+ repmat(pp1, nPP,1).*repmat(1-plw, 1,nTrials); 
V(:,:,2) = repmat(pl2, nPP,1).*repmat(plw, 1,nTrials)+repmat(pp2, nPP,1).*repmat(1-plw, 1,nTrials);
  
end


