function [d] = f_calcnew( d, col, scorebins,fxns)
% Calculate new scores for session 2/3 choice data 
%   Fundamental scores (from PL & PP scores only) are calculated elsewhere,
%   this script is for everything else 
%

% Conflict
d(:, col.PLvcf) = d(:, col.PLcf).*(d(:, col.PL1).*d(:, col.PL2)./100);
d(:, col.PPvcf) = d(:, col.PPcf).*(d(:, col.PP1).*d(:, col.PP2)./100);
d(:, col.PLPPvcf) = d(:, col.PLPPcf).*(d(:, col.PLPP1).*d(:, col.PLPP2)./200);

% keyboard

% Bin scores (according to request)
if isempty(scorebins)==0
    
    % Bin raw scores 
    for i=1:size( scorebins,1)
        ws.binindex= sum(cell2mat( cellfun(@(x, d)d==x, num2cell(scorebins{i})', repmat({d(:, col.PL(1))},  length(scorebins{i}),1), 'UniformOutput',0)'),2 );
        d(find(ws.binindex), col.PLbin1)=i;
        ws.binindex= sum(cell2mat( cellfun(@(x, d)d==x, num2cell(scorebins{i})', repmat({d(:, col.PL(2))},  length(scorebins{i}),1), 'UniformOutput',0)'),2 );
        d(find(ws.binindex), col.PLbin2)=i;
        ws.binindex= sum(cell2mat( cellfun(@(x, d)d==x, num2cell(scorebins{i})', repmat({d(:, col.PP1)},  length(scorebins{i}),1), 'UniformOutput',0)'),2 );
        d(find(ws.binindex), col.PPbin1)=i;
        ws.binindex= sum(cell2mat( cellfun(@(x, d)d==x, num2cell(scorebins{i})', repmat({d(:, col.PP2)},  length(scorebins{i}),1), 'UniformOutput',0)'),2 );
        d(find(ws.binindex), col.PPbin2)=i;
    end
    d(:,col.PLPPbin1)= d(:,col.PLbin1)+ d(:,col.PPbin1);
    d(:,col.PLPPbin2)= d(:,col.PLbin2)+ d(:,col.PPbin2);
    d(:, col.PLbincf) = fxns.fcf( d(:, col.PLbin1),  d(:, col.PLbin2) ) ; 
    d(:, col.PPbincf) = fxns.fcf(d(:, col.PPbin1),  d(:, col.PPbin2) );
    if size(scorebins,1) ==3; d(:, [col.PLbincf col.PPbincf]) =d(:, [col.PLbincf col.PPbincf]) -14; 
    else error('How much to correct binned scores by?'); 
%         Remember that it needs to accomodate PL/PP alone, and PL+PP
    end
    
    
    d(:, col.PLbinvcf) = d(:, col.PLbincf).*(d(:, col.PLbin1).*d(:, col.PLbin2)) ./10;
    d(:, col.PPbinvcf) = d(:, col.PPbincf).*(d(:, col.PPbin1).*d(:, col.PPbin2)) ./10;
    d(:, col.ad1m2.PLbin) = abs(d(:, col.PLbin1) - d(:, col.PLbin2)) ;
    d(:, col.ad1m2.PPbin) = abs(d(:, col.PPbin1) - d(:, col.PPbin2)) ;
    
    % Re-cast binned scores in terms of Cho vs Uncho
    d(d(:, col.Choice)==1, [col.choPPbin  col.choPLbin  col.choPLPPbin]) =  d(d(:, col.Choice)==1, [col.PPbin1     col.PLbin1      col.PLPPbin1]);
    d(d(:, col.Choice)==2, [col.choPPbin  col.choPLbin   col.choPLPPbin]) = d(d(:, col.Choice)==2, [col.PPbin2      col.PLbin2      col.PLPPbin2]);
    d(d(:, col.Choice)==1, [col.unchoPPbin   col.unchoPLbin   col.unchoPLPPbin]) =  d(d(:, col.Choice)==1, [col.PPbin2          col.PLbin2          col.PLPPbin2]);
    d(d(:, col.Choice)==2, [col.unchoPPbin   col.unchoPLbin   col.unchoPLPPbin]) =  d(d(:, col.Choice)==2, [col.PPbin1          col.PLbin2          col.PLPPbin1]);
    d(:, col.marPL) =    d(:, col.choPL)- d(:, col.unchoPL);
    d(:, col.marPP) =    d(:, col.choPP)- d(:, col.unchoPP);
    d(:, col.marPLPP) =    d(:, col.choPLPP)- d(:, col.unchoPLPP);
    d(:, col.PLmarbin) =    d(:, col.choPLbin)- d(:, col.unchoPLbin);
    d(:, col.PPmarbin) =    d(:, col.choPPbin)- d(:, col.unchoPPbin);
    d(:, col.PLPPmarbin) =    d(:, col.choPLPPbin)- d(:, col.unchoPLPPbin); 
end



end

