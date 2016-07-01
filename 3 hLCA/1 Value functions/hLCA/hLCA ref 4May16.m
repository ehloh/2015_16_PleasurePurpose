function [ y ] = hLCA(d, x, fp) % d=data 
% 

% Free parameters 
x.nd = 7000;   % Non-decision time (ms)

x.decay= -0.1;
x.latinhib= -0.5;
x.threshold= 5; 
x.attribeta= 0.4;  % Attribute beta (in attribute weighting)


fp.noisevar= 0.4; 
fp.timebin =10;  % Size of time bin (ms)

% Fixed parameters in set up 
ntr=size(d,1);
ntb =  floor(x.nd/fp.timebin);  % n time bins 
fvect= @(mat)mat(:);
col=fp.col; 

%%  Feature-level evidence computation  



for tr=1:ntr 
%     d(tr, [col.pl1 col.pp1 col.pl2 col.pp2]) 
    tr_fev= [zeros(1,4); nan(ntb-1,4)]; % Attribute-level
    tr_oev= [zeros(1,2); nan(ntb-1,2)]; % Option-level 
 
    for tb=2:ntb
        
        % Feature-level evidence accumulation 
        tr_fev(tb, :) = tr_fev(tb-1,:) ...                  
            + d(tr, [col.pl1 col.pp1 col.pl2 col.pp2]) ...      % Driving response: attribute values 
            + x.decay.*tr_fev(tb-1, :) ...                        % Decay 
            + x.latinhib.*[tr_fev(tb-1,3:4)  tr_fev(tb-1,1:2)] ... % Lateral inhibition (within attribute)
            + fp.noisevar.*randn(1,4);                              % Noise  
        
        % Attribute-weighting 
        weightPL = 1/(1+ exp(x.attribeta* (abs(tr_fev(tb, 1)-tr_fev(tb, 3)) - abs(tr_fev(tb, 2)-tr_fev(tb, 4)))     ) ); 
          
        % Integrate across features
        tr_oev(tb, :) = tr_oev(tb-1,:) ...
            + [sum(tr_fev(tb, [1 2]).*[weightPL 1-weightPL])  sum(tr_fev(tb, [3  4]).*[weightPL 1-weightPL])] ...   % Driving response: overall value (integrated across dimensions)
            + x.decay.*tr_oev(tb-1, :) ...                        % Decay
            + x.latinhib.*tr_oev(tb-1, [2 1]) ...                % Lateral inhibition
            + fp.noisevar.*randn(1,2);                            % Noise
        
        % We good? 
        if sum(tr_oev(tb, :) >  x.threshold)>1; break; end  
    end 
    
%     subplot(1,2,1),  imagesc(tr_fev(1:tb,:)), colorbar, title('Feature-level evidence')
%     subplot(1,2,2),  imagesc(tr_oev(1:tb,:)), colorbar, title('Option-level evidence')
      
end





for o=1:1  % Implement in a matrix
    this = 0;
    if this
        ntr=size(d,3);  % n trials
        
        %
        f_ev=zeros(2,2, ntr, ntb);  % Evidence (Initialize at 0)
        f_contr = eye(2);  % Contrast matrix -  not sure I understand this yet. Laurence subtracts the mean - leave it for now.
        f_change=  [x.decay x.latinhib; x.latinhib x.decay];   % Competition and decay
        
        
        
        disp('the problem here is that i couldn''t construct the matrix multiplication w so many dimensions (pl, pp, tr)')
            
        for tb= 2:ntb
           
            %
            %     keyboard
            %
            %     f_ev(:,:,:,tb) = f_ev(:,:,:,tb-1)
            %
            %
            %     +
            %
            %     f_ev(:,:,:,tb-1)'
            %
            %     *
            %
            %     repmat([x.decay  x.latinhib; x.latinhib x.decay], [1 1 ntr])
            %
        end
    end

end 
   



for o=1:1  % Break it apart from the matrix formulation entirely
    this =0;
    if this
        ntr=size(d,1);
        f_ev=zeros(4, ntr, ntb);   % PL1, PP1, PL2, PP2
        
        disp('problem w this is that i cant figure out intuitively what the matrix operations are doing in order to break them down. if i could figure this out i could probably continue w this approach')
        for tb= 2:ntb 
            f_ev(1:4, :, tb)=f_ev(1:4, :, tb-1)...
                + x.decay*f_ev(1:4, :, tb-1) ... % Decay
                + x.latinhib*f_ev(1:4, :, tb-1) ... % Lateral Inhibition
                + fp.noisevar.*randn(1, 4) ... % Noise
                
            fp.noisevar
            
            d(:,  [col.pl1 col.pp1 col.pl2 col.pp2])
        end
    end
    
     
end 

% Implement via the matrix formulae 
ntr=size(d,1);
f_ev=zeros(ntr,4);   % PL1, PP1, PL2, PP2
f_contr = [1 -0.5; -0.5 1];  % Contrast matrix -  not sure I understand this yet. Laurence subtracts the mean - leave it for now.
f_contr =[1 -0.5; 1 -0.5]; 
f_contr =[1 -1; 1 -1]; 
f_contr =[1  1; 1 1];

  


end

