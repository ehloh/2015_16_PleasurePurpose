function [nLL beh] = hLCA(xx, modinputs)   
%   [Inputs]   xx:    free params   
%                  modinputs{1}:   data (defined in fp.col)
%                  modinputs{2}:   fp, fixed parameters/settings for fit procedure 
%                           maxtime: maximum run time (ms)
%                           timebin: length of 1 evidence-collection timebin (ms)
%                           nsamples: n samples to run to calculate p(Observed Choice) - disabled for now 
%   [Outputs] 
%     nLL: nLL for choice 
%     beh: predicted choice, RT (row=trial)
%     simev: simulated evidence (Trial, Timebin, Option)
%  ------------------------------------------------------------------------

d= modinputs{1};  fp= modinputs{2}; 
 
% Free parameters + apply constraints 
x.decay =    -1./(1+exp(-xx(1))); 
x.latinhib =       -1./(1+exp(-xx(2))); 
x.threshold =        100./(1+exp(-xx(3))); 
x.attribeta =        10./(1+exp(-xx(4))); 
x.noisevar  =xx(5);

% Logistics 
ntr=size(d,1);
ntb =  floor(fp.maxtime/fp.timebin);  % n time bins 
% fvect= @(mat)mat(:);
col=fp.col; 

%%  Evidence computation   

beh=nan(ntr, 2, fp.nsamples); % [2]: Choice, RT

for sm= 1: fp.nsamples
    simev= nan(ntr, ntb,2); 
    
    for tr=1:ntr 
    %     d(tr, [col.pl1 col.pp1 col.pl2 col.pp2]) 
        tr_fev= [zeros(1,4); nan(ntb-1,4)]; % Attribute-level
        tr_oev= [zeros(1,2); nan(ntb-1,2)]; % Option-level 

        for tb=2:ntb

            % Feature-level evidence accumulation 
            tr_fev(tb, :) = tr_fev(tb-1,:) ...                  
                - mean(tr_fev(tb-1,:) ) ...                             % Mean centre
                + d(tr, [col.pl1 col.pp1 col.pl2 col.pp2]) ...      % Driving response: attribute values 
                + x.decay.*tr_fev(tb-1, :) ...                        % Decay 
                + x.latinhib.*[tr_fev(tb-1,3:4)  tr_fev(tb-1,1:2)] ... % Lateral inhibition (within attribute)
                + x.noisevar.*randn(1,4);                              % Noise  

            % Attribute-weighting 
            weightPL = 1/( 1+ exp(      x.attribeta* ( abs(tr_fev(tb, 1)-tr_fev(tb, 3)) - abs(tr_fev(tb, 2)-tr_fev(tb, 4)) )     ) ); 

            % Integrate across features
            tr_oev(tb, :) = tr_oev(tb-1,:) ...
                - mean(tr_oev(tb-1,:) ) ...                             % Mean centre
                + [sum(tr_fev(tb, [1 2]).*[weightPL 1-weightPL])  sum(tr_fev(tb, [3  4]).*[weightPL 1-weightPL])] ...   % Driving response: overall value (integrated across dimensions)
                + x.decay.*tr_oev(tb-1, :) ...                        % Decay
                + x.latinhib.*tr_oev(tb-1, [2 1]) ...                % Lateral inhibition
                + x.noisevar.*randn(1,2);                              % Noise  

            % We good? 
            if any(tr_oev(tb, :) >  x.threshold); break; 
%             elseif tb==ntb,  disp('Reached end of evidence collection without a decision! Now what?'); keyboard 
    %                 ans = tr_oev
    %                 tr_oev(tb-1,:)
            end   
        end 

        % Output results of model simulations 
         simev(tr, 1:tb, :) =  tr_oev(1:tb, :); % Trial, Timebin, Option 
         if length([find(tr_oev(tb, :)>x.threshold) tb*fp.timebin])~=2,   % Patch: If evidence for both options exceeds in 1 timestep, pick the one w the most evidence 
             beh(tr, 1:2,sm) =  [find(tr_oev(tb, :)==max(tr_oev(tb, :)))   tb*fp.timebin];
         else  beh(tr, 1:2,sm) =  [find(tr_oev(tb, :)>x.threshold) tb*fp.timebin]; 
         end   
         
    %     % %  % Wanna see? 
    %     subplot(1,2,1),   plot(1:tb,  tr_fev(1:tb,:),'linewidth',2), legend('PL1','PP1', 'PL2','PP2');  hold on, plot(0:tb,  repmat(x.threshold,1,tb+1), 'color','k'), title('Feature-level evidence', 'FontSize',20)
    %     set(gca, 'xtick', 2:tb, 'xticklabel', (2:tb)*fp.timebin, 'FontSize',20), xlabel('ms'), ylabel('Evidence') 
    %     subplot(1,2,2),   plot(1:tb,  tr_oev(1:tb,:),'linewidth',2), legend('Opt1', 'Opt2');  hold on, plot(0:tb,  repmat(x.threshold,1,tb+1), 'color','k'),  title('Option-level evidence', 'FontSize',20)
    %     set(gca, 'xtick', 2:tb, 'xticklabel', (2:tb)*fp.timebin, 'FontSize',20), xlabel('ms'), ylabel('Evidence') 
    %     disp(['Choice = ' num2str(find(tr_oev(tb, :)>x.threshold)) '  RT = '  num2str(tb*fp.timebin) ' ms'])
    end    
end 

% Calculate likelihood of observed choice (RTs are ignored for now)
pcho= nan(ntr, 2); 
for tr=1:ntr 
    pcho(tr, 1:2)=  [mean(squeeze(beh(tr, 1, :))==1) mean(squeeze(beh(tr, 1, :))==2)]; 
end
nLL = - sum(log( pcho(  d(:, col.ch)) ));  
 
% cl=clock;  disp(['     exiting model (' num2str(cl(4)) ':'  num2str(cl(5)) ')'])
 

% % Model-predicted choice assessed wrt binomial distribution (e.g. Tsetsos 2012 Frontiers)
% k= sum(beh(:,1)==d(:, col.ch));   
% p=  mean(beh(:,1)==d(:, col.ch) );
% nLL =  -log( (p^k)*(  (1-p)^(ntr-k) )  );     

end

