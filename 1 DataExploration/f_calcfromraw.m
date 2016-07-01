function [ d] = f_calcfromraw( d, col)
% Calculate scores from raw data (PL & PP scores only), for session 2/3 choice data
%   Only the most fundamental scores are calculated here  


% Differences
d(:, col.PLPP(1) ) = d(:, col.PL(1)) + d(:, col.PP(1));
d(:, col.PLPP(2) ) = d(:, col.PL(2)) + d(:, col.PP(2));
d(:, col.d1m2.PL)=  d(:, col.PL(1)) - d(:, col.PL(2));
d(:, col.d1m2.PP)=  d(:, col.PP(1)) - d(:, col.PP(2));
d(:, col.d1m2.PLPP)=  d(:, col.PLPP(1)) - d(:, col.PLPP(2));
d(:, col.ad1m2.PL)=  abs( d(:, col.PL(1)) - d(:, col.PL(2)) );
d(:, col.ad1m2.PP)=  abs( d(:, col.PP(1)) - d(:, col.PP(2)) );
d(:, col.ad1m2.PLPP)=  abs( d(:, col.PLPP(1)) - d(:, col.PLPP(2)) );

% What counts as an incompatible-choice  trial? 
d(d(:, col.PL(1))>d(:, col.PL(2)), col.PLbest) =1;
d(d(:, col.PL(1))<d(:, col.PL(2)), col.PLbest) =2;
d(d(:, col.PP(1))>d(:, col.PP(2)), col.PPbest) =1;
d(d(:, col.PP(1))<d(:, col.PP(2)), col.PPbest) =2; 
d(d(:, col.PLPP(1))>d(:, col.PLPP(2)), col.PLPPbest) =1;
d(d(:, col.PLPP(1))<d(:, col.PLPP(2)), col.PLPPbest) =2;
d(:, col.Trialcf)= d(:, col.PLbest)~=d(:, col.PPbest);
d(d(:, col.PLbest)==0 | d(:, col.PPbest)==0, col.Trialcf)=0; % Draws ~=conflict 

% Cho/Uncho
d(d(:, col.Choice)==1, [col.choPP col.choPL col.choPLPP]) = d(d(:, col.Choice)==1,  [col.PP(1) col.PL(1) col.PLPP(1)]);
d(d(:, col.Choice)==1, [col.unchoPP col.unchoPL col.unchoPLPP]) =d(d(:, col.Choice)==1,  [col.PP(2) col.PL(2) col.PLPP(2)]);
d(d(:, col.Choice)==2, [col.choPP col.choPL col.choPLPP]) =d(d(:, col.Choice)==2,  [col.PP(2) col.PL(2) col.PLPP(2)]);
d(d(:, col.Choice)==2, [col.unchoPP col.unchoPL col.unchoPLPP]) =d(d(:, col.Choice)==2,  [col.PP(1) col.PL(1) col.PLPP(1)]);

% Choice congruent with PL/PP/PLPP
%       Draws are completed uncounted (nans)
d(:, [col.conchoPL col.conchoPP col.conchoPLPP])=nan;
d( d(:, col.choPL) >d(:, col.unchoPL) , col.conchoPL) =1;    
d( d(:, col.choPL)<d(:, col.unchoPL) , col.conchoPL) =-1;
d( d(:, col.choPP) >d(:, col.unchoPP) , col.conchoPP) =1;
d( d(:, col.choPP)<d(:, col.unchoPP) , col.conchoPP) =-1;
d( d(:, col.choPLPP) >d(:, col.unchoPLPP) , col.conchoPLPP) =1;
d( d(:, col.choPLPP)<d(:, col.unchoPLPP) , col.conchoPLPP) =-1;


end

