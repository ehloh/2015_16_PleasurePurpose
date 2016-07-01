% Flexibly plot/explore data 
clear all;close all hidden; clc;   logg.subsample_trials=[];  request.session=2;    % 2= fMRI, 3= in lab 
where.fmri_scripts ='/Users/EleanorL/Dropbox/SCRIPPS/3b PLPP fmri';  addpath(where.fmri_scripts)
   
for o=1:1  % Subject subsets 
    % Median split by PL-PP correlation 
    logg.samples.plppcorlow= {'p04_AW';'p12_AL';'p09_KJ';'p01_YH';'p11_BI';'p08_KC';'p17_BB';'p18_WB';'p19_HB';'p02_MI';}; 
    logg.samples.plppcorhigh={'p14_SK';'p15_MM';'p06_BT';'p10_YC';'p13_MS';'p16_SH';'p03_AY';'p07_HC';'p05_CA';'p20_LZ';};
    logg.samples.cc_plppcorlow={'p12_AL';'p17_BB';'p19_HB';'p11_BI';'p08_KC';'p04_AW';'p01_YH';'p18_WB';'p15_MM';'p14_SK';}; 
    logg.samples.cc_plppcorhigh={'p02_MI';'p06_BT';'p09_KJ';'p10_YC';'p13_MS';'p03_AY';'p16_SH';'p07_HC';'p05_CA';'p20_LZ';};
    logg.samples.ic_plppcorlow= {'p19_HB';'p09_KJ';'p17_BB';'p11_BI';'p04_AW';'p01_YH';'p12_AL';'p18_WB';'p08_KC';'p02_MI';}; 
    logg.samples.ic_plppcorhigh= {'p14_SK';'p15_MM';'p06_BT';'p07_HC';'p16_SH';'p05_CA';'p10_YC';'p13_MS';'p20_LZ';'p03_AY';};
    logg.samples.pilot.plppcorlow={'p11';'p04';'p13';'p06';'p19';'p07';'p10';'p12';'p18';'p09';}; 
    logg.samples.pilot.plppcorhigh={'p02';'p05';'p16';'p01';'p15';'p20';'p17';'p03';'p08';'p14';}; 
    logg.samples.pilot.cc_plppcorlow={'p13';'p11';'p06';'p04';'p07';'p10';'p09';'p12';'p02';'p16';};
    logg.samples.pilot.cc_plppcorhigh={'p18';'p19';'p05';'p01';'p15';'p17';'p03';'p20';'p08';'p14';}; 
    logg.samples.pilot.ic_plppcorlow= {'p18';'p11';'p19';'p04';'p06';'p07';'p10';'p09';'p13';'p12';}; 
    logg.samples.pilot.ic_plppcorhigh={'p05';'p02';'p01';'p16';'p17';'p14';'p20';'p08';'p15';'p03';}; 
    %
    logg.samples.nd_plppcorlow= {'p19_HB';'p08_KC';'p01_YH';'p04_AW';'p17_BB';'p12_AL';'p11_BI';'p09_KJ';'p02_MI';'p18_WB';};
    logg.samples.nd_plppcorhigh= {'p14_SK';'p15_MM';'p06_BT';'p07_HC';'p16_SH';'p03_AY';'p10_YC';'p05_CA';'p13_MS';'p20_LZ';};
    logg.samples.nd_cc_plppcorlow= {'p19_HB';'p08_KC';'p01_YH';'p04_AW';'p17_BB';'p12_AL';'p11_BI';'p09_KJ';'p02_MI';'p18_WB';};
    logg.samples.nd_cc_plppcorhigh= {'p14_SK';'p15_MM';'p06_BT';'p07_HC';'p16_SH';'p03_AY';'p10_YC';'p05_CA';'p13_MS';'p20_LZ';};
    %
    logg.samples.ppresid_plppcorlow= {'p09_KJ';'p01_YH';'p04_AW';'p12_AL';'p02_MI';'p11_BI';'p08_KC';'p18_WB';'p19_HB';'p14_SK';};
    logg.samples.ppresid_plppcorhigh= {'p10_YC';'p17_BB';'p15_MM';'p07_HC';'p06_BT';'p16_SH';'p13_MS';'p05_CA';'p20_LZ';'p03_AY';};
    logg.samples.ppresid_cc_plppcorlow={'p12_AL';'p11_BI';'p10_YC';'p14_SK';'p19_HB';'p02_MI';'p04_AW';'p15_MM';'p16_SH';'p18_WB';};
    logg.samples.ppresid_cc_plppcorhigh={'p01_YH';'p20_LZ';'p07_HC';'p03_AY';'p08_KC';'p09_KJ';'p06_BT';'p13_MS';'p05_CA';'p17_BB';};
    logg.samples.ppresid_ic_plppcorlow={'p09_KJ';'p01_YH';'p02_MI';'p08_KC';'p18_WB';'p19_HB';'p11_BI';'p04_AW';'p17_BB';'p14_SK';};
    logg.samples.ppresid_ic_plppcorhigh={'p12_AL';'p07_HC';'p10_YC';'p06_BT';'p15_MM';'p05_CA';'p16_SH';'p20_LZ';'p13_MS';'p03_AY';};
    
    
    logg.samples.not_p19={'p01_YH';'p02_MI';'p03_AY';'p04_AW';'p05_CA';'p06_BT';'p07_HC';'p08_KC';'p09_KJ';'p10_YC';'p11_BI';'p12_AL';'p13_MS';'p14_SK';'p15_MM';'p16_SH';'p17_BB';'p18_WB';'p20_LZ';};
    % logg.samples.pPLunder90={'p01_YH';'p02_MI';'p03_AY';'p05_CA';'p06_BT';'p07_HC';'p09_KJ';'p10_YC';'p11_BI';'p12_AL';'p13_MS';'p14_SK';'p16_SH';'p17_BB';'p18_WB';'p20_LZ';};   % Exclude people who choose higher PL scores >90% of the time 
    
    
end

% Request specific
for o=1:1  % Actively request subsampels 
    logg.specificsubjects=unique([
        logg.samples.plppcorlow ;
        logg.samples.plppcorhigh ;
        logg.samples.cc_plppcorlow;
        logg.samples.cc_plppcorhigh ;
        logg.samples.ic_plppcorlow;
        logg.samples.ic_plppcorhigh;
        logg.samples.pilot.plppcorlow;
        logg.samples.pilot.plppcorhigh;
        logg.samples.pilot.cc_plppcorlow;
        logg.samples.pilot.cc_plppcorhigh;
        logg.samples.pilot.ic_plppcorlow;
        logg.samples.pilot.ic_plppcorhigh;
        ]);
end
logg.specificsubjects={};    
% logg.specificsubjects={'p01_YH'};  
% logg.specificsubjects={'p10_YC'};


% logg.specificsubjects= f_subsample('beh1','PLvPP_abs_r04'); 


% Dataset 
% request.prebinned =' b3'; 
request.prebinned =[]; 
request.dataset=['All data' request.prebinned ' (17-Nov-2015)'];   request.session=2;  
% request.dataset=['All pilot data' request.prebinned ' (17-Nov-2015)'];  request.session=1;  



% request.dataset=['All data' request.prebinned ' (26-May-2016)'];   request.session=2;    % Different ordinal bins (middle is 4)

for o1=1:1 % General settings and specifications
   
    % Load subjects
    where.where='/Users/EleanorL/Dropbox/SCRIPPS/3 Pleasure purpose';  
    where.beh=[where.where filesep '3 Behaviour'];    where.mod =[where.where filesep '4b Modelling']; 
    path(pathdef), addpath(where.where), addpath(where.mod)  
    addpath([where.mod filesep '1 Value functions'])  
    addpath([where.mod filesep '1 Value functions' filesep 'hLCA'])  
    if ischar(logg.specificsubjects)==1 && strcmp(logg.specificsubjects, 'AllDataOK'), logg.specificsubjects= {'p01_YH';'p02_MI';'p03_AY';'p04_AW';'p05_CA';'p06_BT';'p07_HC';'p08_KC';'p09_KJ';'p10_YC';'p11_BI';'p13_MS';'p14_SK';'p15_MM';'p16_SH';'p17_BB';'p18_WB';'p19_HB';'p20_LZ';}; end ;
    if isempty(strfind(request.dataset, 'pilot'))==1 & isempty(strfind(request.dataset, 'Pilot'))==1
        [n t r]=xlsread([where.beh filesep 'datalog_plpr.xlsx']);
        [logg.subjects logg.n_subjs logg.datalog] = f_selectsubjects(r, logg.specificsubjects,r, 'All');
    else
        logg.subjects=[cellfun(@(x)['p0' num2str(x)], num2cell(1:9),'UniformOutput',0)'; cellfun(@(x)['p' num2str(x)], num2cell(10:20),'UniformOutput',0)';];
        logg.n_subjs = length(logg.subjects); 
        [logg.subjects logg.n_subjs logg.datalog] = f_selectsubjects([[{'Subject'}; logg.subjects]  [{'All'}; num2cell(ones(20,1))] ], logg.specificsubjects,[[{'Subject'}; logg.subjects]  [{'All'}; num2cell(ones(20,1))] ], 'All');
    end
    w=load([where.where fs '4a Beh analysis basic' fs '2 Data' fs request.dataset '.mat']); 
    subjdata=w.subjdata; col=w.col; logg.define =  w.log.define; logg.fxns.fcf= logg.define{strcmp(logg.define(:,1), 'Conflict'),3}; 
    [d d subjdata] = f_selectsubjects([{'Subject'}  cell(1, size(subjdata,2)-1 )  ; subjdata ], logg.subjects,[[{'Subject'}  cell(1, size(subjdata,2)-1 ) {'All'} ]; [subjdata num2cell(ones(size(subjdata,1),1))]], 'All');   subjdata=subjdata(2:end, :);
   
    % Misc
    request.alliv_names= {'PL1' 'PL1'; 'PL2' 'PL2';'PP1' 'PP1'; 'PP2' 'PP2';  'd1m2.PL' 'PL diff'; 'd1m2.PP' 'PP diff';    'PLcf' 'PL cf';  'PPcf' 'PP cf';   'PLPP1' 'PLPP1';   'PLPP2' 'PLPP2';  'd1m2.PLPP' 'PLPP diff';       'PLPPcf' 'PLPP cf'};       
    logg.define={}; 
    eval(['col=col.s' num2str(request.session) ';'])  
    logg.define(size(logg.define,1)+1, :)={'RT preproc' 'Outliners (>2SD +/- mean) are nan-ed' [] };
    request.IVs_rescale =0;  % dont

    % Interface=
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
    disp(['No. of subjects: ' num2str(logg.n_subjs)])
    if isempty(logg.specificsubjects)==0;  if logg.n_subjs == 19 & sum(strcmp(logg.subjects, {'p01_YH';'p02_MI';'p03_AY';'p04_AW';'p05_CA';'p06_BT';'p07_HC';'p08_KC';'p09_KJ';'p10_YC';'p11_BI';'p13_MS';'p14_SK';'p15_MM';'p16_SH';'p17_BB';'p18_WB';'p19_HB';'p20_LZ';}) ) ==19; disp('   Subjects w complete data only (n=19)'); else disp('   Subset of subjects only:'); disp(logg.specificsubjects); end;   end; disp(' ')
    disp(' '); input('Hit Enter to start      ')
    disp('=======================================================')
    
end

%%  Redo choice via simulation from the hLCA


request.do_sim=1;   
mod.name='lta01';
if request.do_sim
    for o=1:1  % Settings for simulation 
        mod.req.normalize_withinattr =1;
        mod.fp.nsamples=100;
        %
        mod.fp.maxtime= 7000;
        mod.fp.timebin =10;  % Size of time bin (ms)
        
        % Settings for free parameters
        [mod.defaults  mod.par mod.models] = f_hLCAsettings({mod.name}, 10);
        mod.setup =  [mod.models{3}' mod.models{4}' mod.models{5}' mod.models{6}' mod.models{7}'];
        for p= 1: length(mod.defaults{ strcmp(mod.defaults(:,1), 'hLCA'),3}) % Assign defaults for key params
            if ~isempty(mod.par{strcmp(mod.par(:,1), mod.defaults{ strcmp(mod.defaults(:,1), 'hLCA'),3}{p}),5} ),  eval(['mod.fp.'  mod.defaults{ strcmp(mod.defaults(:,1), 'hLCA'),3}{p} ' = '  num2str(mod.par{strcmp(mod.par(:,1), mod.defaults{ strcmp(mod.defaults(:,1), 'hLCA'),3}{p}),5}) ';']), end
        end
        
        % Columns for fed-in data
        mod.fp.col.pl1=1;
        mod.fp.col.pp1=2;
        mod.fp.col.pl2=3;
        mod.fp.col.pp2=4;
        mod.fp.col.ch=5;
        mod.fp.col.rt=6;
        
    end
    for s=1:logg.n_subjs
        disp([logg.subjects{s} '  -  Simulating choices'])
        ws.dd= subjdata{s, request.session+1};
        
        % Get rid of trials in which there are no ratings
        ws.dd = ws.dd(~isnan(sum(ws.dd(:, [ col.PL(1) col.PP(1) col.PL(2) col.PP(2)]),2)),:);
        subjdata{s, request.session+1} = ws.dd;
        
        if mod.req.normalize_withinattr % Normalize within-attribute
            ws.dd(:, [ col.PL(1) col.PL(2)]) = ws.dd(:, [ col.PL(1) col.PL(2)])./repmat(sum(ws.dd(:, [ col.PL(1) col.PL(2)]),2),1,2);
            ws.dd(:, [ col.PP(1) col.PP(2)]) = ws.dd(:, [ col.PP(1) col.PP(2)])./repmat(sum(ws.dd(:, [ col.PP(1) col.PP(2)]),2),1,2);
        end
        ws.d= ws.dd(:, [col.PL(1) col.PP(1) col.PL(2) col.PP(2) col.Choice col.ChoiceRT]);
        ws.d(:, mod.fp.col.ch)=randi(2,size(ws.d,1),1);
        
        % Set up parameters
        for p=1: size(mod.par,1)
            x  = mod.par{p,4}(round(length(mod.par{p,4})/2));  % Midpoint  of range
            eval(['ws.par(p)='  mod.par{p,3} ';'])
        end
        eval(['[ws.nLL ws.beh] =' mod.name  '(ws.par, {ws.d mod.fp}); ']) 
        d_modsim(s,1)=ws.nLL; 

        % Export likely choice
        ws.pcho = nan(size(ws.d,1), 2);
        for t=1:size(ws.d,1)
            ws.pcho(t, 1:2)=   [mean(squeeze(ws.beh(t, 1, :))==1) mean(squeeze(ws.beh(t, 1, :))==2)] ;
            
            if ws.pcho(t, 1)==ws.pcho(t, 2), ws.pcho(t,3)=randi(2); % Draw = random choice
            else  ws.pcho(t,3)=  find(ws.pcho(t, 1:2)==max(ws.pcho(t, 1:2)));
            end
            
            % RT
            ws.pcho(t,4)=  mean( squeeze(ws.beh(t, 2, :) ) );
        end
        
        subjdata{s, request.session+1}(:, col.Choice)= ws.pcho(:,3);
        subjdata{s, request.session+1}(:, col.ChoiceRT)= ws.pcho(:,4);
        
        
        
        
%         
%         
        %% Hard-coded simulations  
        
%         % [Choose 1 dimension randomly and pick the best ] ##############################
%         ws.dd=  ws.dd(ws.dd(:, col.PL(1))~=ws.dd(:, col.PL(2)),:); ws.dd=  ws.dd(ws.dd(:, col.PP(1))~=ws.dd(:, col.PP(2)),:);  if s==1, disp('SIM: Draws deleted in simulation'); end          
%         ws.dd=  repmat(ws.dd, 20,1);  if s==1, disp('SIM: Upsampling '); end           
%         ws.dim = randi(2, size(ws.dd,1),1);  % 1= PL, 2= PP   
%         ws.d1= ws.dd(ws.dim==1, :); ws.d2= ws.dd(ws.dim==2, :);   
%         if s==1, f_sigma=inline('1./(1+exp(-x))');ws.beta= 0.3;  end; 
%         ws.d1p1 =f_sigma(ws.beta*(ws.d1(:,  col.PL(1)) - ws.d1(:,  col.PL(2))));  % Stochasitcally choose best 
%         ws.d1dice = rand(size(ws.d1p1 ,1),1);
%         ws.d1(:, col.Choice) = 2;  ws.d1(ws.d1p1 >ws.d1dice, col.Choice) = 1;  
%         ws.d2p1 =f_sigma(ws.beta*(ws.d2(:,  col.PP(1)) - ws.d2(:,  col.PP(2))));
%         ws.d2dice = rand(size(ws.d2p1 ,1),1);
%         ws.d2(:, col.Choice) = 2;  ws.d2(ws.d2p1 >ws.d2dice, col.Choice) = 1;  
%         ws.dd= [ws.d1; ws.d2]; 
%        
%         
%         
% %          % [Pick PL w noise] ##############################
%         ws.dd=  ws.dd(ws.dd(:, col.PL(1))~=ws.dd(:, col.PL(2)),:); ws.dd=  ws.dd(ws.dd(:, col.PP(1))~=ws.dd(:, col.PP(2)),:);  if s==1, disp('SIM: Draws deleted in simulation'); end          
%         ws.dd=  repmat(ws.dd, 20,1);  if s==1, disp('SIM: Upsampling '); end            
%         ws.d1= ws.dd; 
%         if s==1, f_sigma=inline('1./(1+exp(-x))');ws.beta= 0.3;  end;
%         ws.d1p1 =f_sigma(ws.beta*(ws.d1(:,  col.PL(1)) - ws.d1(:,  col.PL(2))));  % Stochasitcally choose best 
%         ws.d1dice = rand(size(ws.d1p1 ,1),1);
%         ws.d1(:, col.Choice) = 2;  ws.d1(ws.d1p1 >ws.d1dice, col.Choice) = 1;             
%         ws.dd= ws.d1;  
        
        %
        subjdata{s, request.session+1}  =ws.dd;
        subjdata{s, request.session+1}  = f_calcfromraw(subjdata{s, request.session+1}, col);
    end
end
  
%% Calculate new scores 

for o=1:1 % Archive of un-used 
    f_rescale=@(rmin, rmax,x)    rmin + (rmax-rmin).*(x- min(x(:))) ./ max(x- min(x(:))); 
    
    % Note: conflict is already computed by this point. If you want to play w other formulations, 
    % give them other names. Note also that if you want to set up other IVs that are 
    % not choice invariant (i.e. dependent on chosen option's value, not just available options' values),
    % the RT analysis needs to be flipped.
    fsigma=inline('1./(1+exp(-x))');
end
frt=@(x) log(x); logg.define(size(logg.define,1)+1, :)={'RT preproc' 'Log-transformed within subject' frt};
% frt=@(x) zscore(x);   logg.define(size(logg.define,1)+1, :)={'RT preproc' 'Zscored within subject' frt}; disp('RT z scored rather than logged!!!');
    
%     logg.scorebins={[1 2]; [3 4]; [5 6]; [7 8]; [9 10]}; 
    logg.scorebins={[1 2 3]; [4 5  6  7];  [8 9 10]}; 

% Calculate new scores  
for o=1:1 % Columns 
    col.PL1=col.PL(1); col.PL2=col.PL(2);  col.PP1=col.PP(1); col.PP2=col.PP(2); 
    
    col.colmax=structmax(col);    % (generally for IVs, put PL/PP as prefix - for cols)
    col.Subject=col.colmax+1;
    col.PLvcf=col.colmax+5;    
    col.PPvcf=col.colmax+6;    
    col.PLPPvcf=col.colmax+7;    
    col.PLbin1=col.colmax+8;          % Bins of PL/PP scores 
    col.PLbin2=col.colmax+9;    
    col.PPbin1=col.colmax+10;
    col.PPbin2=col.colmax+11;    
    col.PLPPbin1=col.colmax+12;
    col.PLPPbin2=col.colmax+13;     
    %
    col.PLbincf=col.colmax+14;   
    col.PPbincf=col.colmax+15;   
    col.PLPPbincf=col.colmax+16;   
    col.PLbinvcf=col.colmax+17;   
    col.PPbinvcf=col.colmax+18; 
    col.PLPPbinvcf=col.colmax+19;  
    col.ad1m2.PLbin=col.colmax+20;   
    col.ad1m2.PPbin=col.colmax+21; 
    col.ad1m2.PLPPbin=col.colmax+22; 
    col.PLdraw=col.colmax+23; 
    col.PPdraw=col.colmax+24; 
    col.PLPPdraw=col.colmax+25; 
    col.cfweight=col.colmax+26;      % Ratio of PL:PP cf 
    col.cfweightbin=col.colmax+27; 
    col.cfweightlog=col.colmax+28; 
    col.cfweightbinrecip=col.colmax+29; 
    col.nTrials=col.colmax+30; 
    
    col.colmax2=structmax(col);     % Cho, Uncho 
    col.choPLbin= col.colmax2+1;     
    col.choPPbin= col.colmax2+2;     
    col.choPLPPbin= col.colmax2+3;     
    col.unchoPLbin= col.colmax2+4;     
    col.unchoPPbin= col.colmax2+5;     
    col.unchoPLPPbin= col.colmax2+6;     
    col.marPL= col.colmax2+7;  col.PLmar=col.marPL; 
    col.marPP= col.colmax2+8;  col.PPmar=col.marPP;  
    col.marPLPP= col.colmax2+9; 
    col.PLmarbin= col.colmax2+10; 
    col.PPmarbin= col.colmax2+11;  
    col.PLPPmarbin= col.colmax2+12;  
    col.chobest=col.colmax2+13;   % Chose best, assuming equal weighting of attributes? 
    col.pChoPLbin=col.colmax2+14;  
    col.pChoPPbin=col.colmax2+15;  
    col.pChoPLPPbin=col.colmax2+16;  
    
    
    col.colmax3=structmax(col);     % Behaviour (generally, put PL/PP as suffix - for cols)
    col.pChoPL= col.colmax3+1;                      % pCho is per IV bin
    col.pChoPP= col.colmax3+2;   
    col.pChoPLPP= col.colmax3+3;    
    col.RTchoPL= col.colmax3+4;     
    col.RTchoPP= col.colmax3+5;  
    col.RTchoPLPP= col.colmax3+6;  
     
    col.PLad1m2=col.ad1m2.PL;  col.PPad1m2=col.ad1m2.PP;  col.PLbinad1m2=col.ad1m2.PLbin;  col.PPbinad1m2=col.ad1m2.PPbin; 
    col.PLd1m2=col.d1m2.PL;  col.PPd1m2=col.d1m2.PP;   % col.PLbind1m2=col.d1m2.PLbin;  col.PPbind1m2=col.d1m2.PPbin; 
end

subjdata{logg.n_subjs+1,request.session+1}=[];
d_ic=subjdata;d_cc=d_ic; d_cc_nd=d_ic;  d_cc_drawPL=d_cc; d_cc_drawPP=d_cc; d_nd=d_ic; d_dr=d_ic;  d_nplppd=d_ic;  d_nplppd_cc=d_ic; d_nplppd_ic=d_ic;
for s=1:logg.n_subjs
    ws.d= subjdata{s,request.session+1}; ws.d(:, col.Subject)=s;  ws.d = ws.d(ws.d(:, col.TrialOK)==1,:);
    
%     if s==1; disp('Artificial trial order'),  ws.d=sortrows(ws.d,  [col.PL1 col.PP1 col.PP2 col.PL2 ]); end
%         disp('Hardcode again cf trials!'), ws.d((ws.d(:, col.PL(1))>ws.d(:, col.PL(2))) & (ws.d(:, col.PP(1))<ws.d(:, col.PP(2))), col.Trialcf)=1;  ws.d((ws.d(:, col.PL(1))<ws.d(:, col.PL(2))) & (ws.d(:, col.PP(1))>ws.d(:, col.PP(2))), col.Trialcf)=1;             

    % PP = residuals 
%     if s==1,  input('PP = has had PL regressed out (rescal+rounded). Continue? ');  logg.define(size(logg.define,1)+1,1:3)= {'PP has had PL regressed out (rescaled 1-10 and rounded)' ' ' ' '};  end;     [ws.b ws.bi ws.r]=  regress(ws.d(:, col.PP1),ws.d(:, col.PL1)); ws.d(:, col.PP1)= round( f_rescale(1, 10, ws.r)); [ws.b ws.bi ws.r]=  regress(ws.d(:, col.PP2),ws.d(:, col.PL2)); ws.d(:, col.PP2)=  round(f_rescale(1, 10, ws.r)); [ ws.d] = f_calcfromraw( ws.d, col) ;    
    ws.d(:, col.PLdraw)= ws.d(:, col.PL(1)) ==ws.d(:, col.PL(2)); 
    ws.d(:, col.PPdraw)= ws.d(:, col.PP(1)) ==ws.d(:, col.PP(2));     

    % % [Deterministic artificial data]
%     if s==1,  input('Fake random data! Continue? '); end; ws.d(:, col.Choice)= randi(2, size(ws.d,1),1);[ ws.d] = f_calcfromraw( ws.d, col) ;
%     if s==1,  input('[FAKE DATA] Choose better PL.  Continue? '); end;  ws.d(:, col.Choice)=  randi(2, size(ws.d,1),1);  ws.d(ws.d(:, col.PL1)>ws.d(:, col.PL2), col.Choice)=1;  ws.d(ws.d(:, col.PL1)<ws.d(:, col.PL2), col.Choice)=2;    ws.d(ws.d(:, col.PL1)==ws.d(:, col.PL2), :)=[]; [ ws.d] = f_calcfromraw( ws.d, col) ;% Delete PL draws
%     if s==1,  input('[FAKE DATA] Choose better PLPP.  Continue? '); end; ws.d(:, col.Choice)=  randi(2, size(ws.d,1),1); ws.d(ws.d(:, col.PLPP1)>ws.d(:, col.PLPP2), col.Choice)=1;  ws.d(ws.d(:, col.PLPP1)<ws.d(:, col.PLPP2), col.Choice)=2; [ ws.d] = f_calcfromraw( ws.d, col) ;
%     if s==1,  input('[FAKE DATA] Choose better cf-weighted PLPP. Continue? '); end;  ws.plpp1_wcf= (ws.d(:, col.PL1)./ws.d(:, col.PLcf))+ (ws.d(:, col.PP1)./ ws.d(:, col.PPcf)); ws.plpp2_wcf=  (ws.d(:, col.PL2)./ws.d(:, col.PLcf))+ (ws.d(:, col.PP2)./ ws.d(:, col.PPcf));  ws.d(:, col.Choice)=  randi(2, size(ws.d,1),1);  ws.d(ws.plpp1_wcf>ws.plpp2_wcf, col.Choice)=1;   ws.d(ws.plpp1_wcf<ws.plpp2_wcf, col.Choice)=2; [ ws.d] = f_calcfromraw( ws.d, col) ;

    % % [Stochastic artificial data]
%     ws.beta=0.8*1; ws.d=repmat(ws.d,1000,1);  ws.pl= ws.d(:, [col.PL1 col.PL2]); ws.pp= ws.d(:, [col.PP1 col.PP2]); ws.cf= ws.d(:, [col.PLcf col.PPcf]);
%     if s==1,  input('[FAKE DATA] Choose better PL.  Continue? '); end; ws.scores= ws.pl;
% %     if s==1,  input('[FAKE DATA] Choose better PLPP.  Continue? '); end; ws.scores= ws.pl+ ws.pp;
% %     if s==1,  input('[FAKE DATA] Choose better cf-weighted PLPP.  Continue? '); end; ws.scores= ws.pl./ws.cf(1)+ ws.pp./ws.cf(2);
%     ws.pch= exp(ws.beta*ws.scores)./repmat(sum(exp(ws.beta*ws.scores),2),1,2);      if sum(isnan(ws.pch(:)))+sum(isinf(ws.pch(:)))>0;   ws.d(isnan(ws.pch(:,1))+ isnan(ws.pch(:,2))==2,:)=[];  ws.scores(isnan(ws.pch(:,1))+ isnan(ws.pch(:,2))==2,:)=[];  ws.pch(isnan(ws.pch(:,1))+ isnan(ws.pch(:,2))==2,:)=[];      if sum(ws.scores(isnan(ws.pch(:,1)),1)<=ws.scores(isnan(ws.pch(:,1)),2)) + sum(ws.scores(isnan(ws.pch(:,2)),2)<=ws.scores(isnan(ws.pch(:,2)),1))>0;  disp('simulation pch nan where it shouldnt default to 1. ugh no idea why. '); end ; ws.pch=ws.pch(:);  ws.pch(isnan(ws.pch))=1;  ws.pch=reshape(ws.pch, size(ws.scores)); end % Assume nans are because the value is too high. 
%     ws.d(:, col.Choice)=2;   ws.d(rand(size(ws.d,1),1)<ws.pch(:,1), col.Choice)=1;     disp(['Subject ' num2str(s)] );  [ ws.d] = f_calcfromraw( ws.d, col) ;    % If using fake data, re-calculate scores that are dependent on choice
%  


    
%     fcf= @(x,y)(5 - abs(x- y) );   logg.fxns.fcf=fcf;  % ONLY for bincf 

    % ------------------------------------------------------------------------------------------------------
    
    % Calculate new scores 
	[ws.d] = f_calcnew(ws.d, col, logg.scorebins, logg.fxns);
    
    % RT pre-processing
    ws.rtout_high = (ws.d(:, col.ChoiceRT) > mean(ws.d(:, col.ChoiceRT)) + 2*std(ws.d(:, col.ChoiceRT))   )   ;
    ws.rtout_low = (ws.d(:, col.ChoiceRT) < mean(ws.d(:, col.ChoiceRT)) - 2*std(ws.d(:, col.ChoiceRT)) )   ;
    ws.d( (ws.rtout_high+ws.rtout_low)>0 , col.ChoiceRT)  =nan;
    ws.d(find(1- (ws.rtout_high+ws.rtout_low)), col.ChoiceRT) = frt(ws.d(find(1- (ws.rtout_high+ws.rtout_low)), col.ChoiceRT) );   % Preproc RTs 
    ws.d(:, [col.RTchoPL col.RTchoPP col.RTchoPLPP])=nan;
    ws.d(ws.d(:, col.Trialcf)==1 & ws.d(:, col.conchoPL)==1,  col.RTchoPL)= ws.d(ws.d(:, col.Trialcf)==1 & ws.d(:, col.conchoPL)==1,  col.ChoiceRT);
    ws.d(ws.d(:, col.Trialcf)==1 & ws.d(:, col.conchoPP)==1,  col.RTchoPP)= ws.d(ws.d(:, col.Trialcf)==1 & ws.d(:, col.conchoPP)==1,  col.ChoiceRT);
    ws.d(ws.d(:, col.Trialcf)==1 & ws.d(:, col.conchoPLPP)==1,  col.RTchoPLPP)= ws.d(ws.d(:, col.Trialcf)==1 & ws.d(:, col.conchoPLPP)==1,  col.ChoiceRT);


    
%     ws.d = sortrows(ws.d, [col.Trialcf col.PLbest col.PPbest]);
%     ws.d(:, [col.Trialcf col.PLbest col.PPbest])

    % ------------------------------------------------------------------------------------------------------
    
    % Output 
    subjdata{logg.n_subjs+1,request.session+1}=[subjdata{logg.n_subjs+1,request.session+1}; ws.d];    
    subjdata{s,request.session+1}=ws.d;
    
    % Special samples
    d_cc{s,request.session+1} = ws.d(ws.d(:, col.Trialcf)==0,:);
    d_cc{logg.n_subjs+1,request.session+1}=[d_cc{logg.n_subjs+1,request.session+1};  d_cc{s,request.session+1}  ];    
    %
    ws.dncf=ws.d(ws.d(:, col.Trialcf)==0,:);
    d_cc_nd{s,request.session+1}= ws.dncf(ws.dncf(:, col.PLdraw)==0 & ws.dncf(:, col.PPdraw)==0, :);
    %     d_cc_nd{s,request.session+1} = ws.d(ws.d(:, col.Trialcf)==0 & ws.d(:, col.PLbest)~=0 & ws.d(:, col.PPbest)~=0 ,:);
    d_cc_nd{logg.n_subjs+1,request.session+1}=[d_cc_nd{logg.n_subjs+1,request.session+1};  d_cc_nd{s,request.session+1}  ];    
    %
    d_ic{s,request.session+1} = ws.d(ws.d(:, col.Trialcf)==1,:);
    d_ic{logg.n_subjs+1,request.session+1}=[d_ic{logg.n_subjs+1,request.session+1};  d_ic{s,request.session+1}  ];    
    %
    d_cc_drawPL{s,request.session+1}  =ws.d(ws.d(:, col.PLbest)==0, :); 
    d_cc_drawPL{logg.n_subjs+1,request.session+1}=[d_cc_drawPL{logg.n_subjs+1,request.session+1};  d_cc_drawPL{s,request.session+1}  ];    
    %
    d_cc_drawPP{s,request.session+1}  =ws.d(ws.d(:, col.PPbest)==0, :); 
    d_cc_drawPP{logg.n_subjs+1,request.session+1}=[d_cc_drawPP{logg.n_subjs+1,request.session+1};  d_cc_drawPP{s,request.session+1}  ];    
    %
    d_nd{s,request.session+1} = ws.d(ws.d(:,  col.PLdraw)==0 & ws.d(:,  col.PPdraw)==0,:);
    d_nd{logg.n_subjs+1,request.session+1}=[d_nd{logg.n_subjs+1,request.session+1};  d_nd{s,request.session+1}  ];    
    %
    d_dr{s,request.session+1}= ws.d(ws.d(:, col.PLdraw)==1 | ws.d(:, col.PPdraw)==1 ,:);
    d_dr{logg.n_subjs+1,request.session+1}=  [d_dr{logg.n_subjs+1,request.session+1};  d_dr{s,request.session+1}];
    %
    d_nplppd{s,request.session+1}= ws.d; 
    d_nplppd{s,request.session+1}(  d_nplppd{s,request.session+1}(:, col.PL1)==d_nplppd{s,request.session+1}(:, col.PP1),:)=[];
    d_nplppd{s,request.session+1}(  d_nplppd{s,request.session+1}(:, col.PL2)==d_nplppd{s,request.session+1}(:, col.PP2),:)=[];
    d_nplppd{logg.n_subjs+1,request.session+1}=[d_nplppd{logg.n_subjs+1,request.session+1};  d_nplppd{s,request.session+1}];
    %
    d_nplppd_cc{s,request.session+1}=d_nplppd{s,request.session+1}(d_nplppd{s,request.session+1}(:, col.Trialcf)==0,:); 
    d_nplppd_ic{s,request.session+1}=d_nplppd{s,request.session+1}(d_nplppd{s,request.session+1}(:, col.Trialcf)==1,:); 
    d_nplppd_ic{logg.n_subjs+1,request.session+1}=[d_nplppd_ic{logg.n_subjs+1,request.session+1};  d_nplppd_ic{s,request.session+1}];
    d_nplppd_cc{logg.n_subjs+1,request.session+1}= [d_nplppd_cc{logg.n_subjs+1,request.session+1};  d_nplppd_cc{s,request.session+1}]; 
      
    ws=[]; 
end; os=subjdata;
  
 
% % [Request: Special samples ONLY?] ##########
% subjdata= d_cc; input('[Subsample] CC (no conflict) trials only (draws not deleted). Continue?  '); logg.subsample_trials='cc';
subjdata= d_ic; input('[Subsample] IC (conflict) trials (intrinsically no draws) only. Continue? ');  logg.subsample_trials='ic';
% subjdata= d_cc_nd;   input('[Subsample] CC trials (NO draws) only. Continue?  ');   logg.subsample_trials='cc_nd';
% subjdata= d_cc_drawPL;   input('[Subsample] CC PL-draw trials only. Continue?  ');   logg.subsample_trials='cc_drawPL';
% subjdata= d_cc_drawPP;   input('[Subsample] CC PP-draw trials only. Continue?  ');  logg.subsample_trials='cc_drawPP';
% subjdata= d_dr;   input('[Subsample] PL or PP Draws only. Continue? ');  logg.subsample_trials='draws';
% subjdata= d_nd;   input('Draws in PL/PP (but not PLPP) deleted. Continue?  ');   logg.subsample_trials='nodraws';
% subjdata= d_nplppd;   input('Events where PL==PP are deleted. Continue?  ');   logg.subsample_trials='delete PL==PP';
% subjdata= d_nplppd_cc;   input('CC trials where PL==PP are deleted. Continue?  ');  logg.subsample_trials='cc, delete PL==PP';
% subjdata= d_nplppd_ic;   input('IC trials where PL==PP are deleted. Continue?  ');  logg.subsample_trials='ic, delete PL==PP';
 % mean(cellfun(@(x)size(x,1), subjdata(1:logg.n_subjs, request.session+1)))    % N trials 
 
for o=1:1 % % % Certain subjects  #######
% req.sample= subjdata(1:end-1,1); % ALL SUBJECTS 
% req.sample= {'p02_MI';    'p17_BB'};  % overall <0.4
% req.sample= [12 17];  % overall <0.4
% req.sample= logg.samples.plppcorlow ;
% req.sample= logg.samples.plppcorhigh ;
% req.sample= logg.samples.cc_plppcorlow;
% req.sample= logg.samples.cc_plppcorhigh ;
% req.sample= logg.samples.ic_plppcorlow;
% req.sample= logg.samples.ic_plppcorhigh; 
% req.sample= logg.samples.nd_plppcorlow ;
% req.sample= logg.samples.nd_plppcorhigh; 
% req.sample= logg.samples.nd_cc_plppcorlow ;
% req.sample= logg.samples.nd_cc_plppcorhigh; 
% req.sample= logg.samples.ppresid_plppcorlow;
% req.sample= logg.samples.ppresid_plppcorhigh;
% req.sample= logg.samples.ppresid_cc_plppcorlow;
% req.sample= logg.samples.ppresid_cc_plppcorhigh;
% req.sample= logg.samples.ppresid_ic_plppcorlow;
% req.sample= logg.samples.ppresid_ic_plppcorhigh;

% whichsubs=sortrows(  cell2mat(cellfun(@(x,y)find(strcmp(x,y)),  repmat({subjdata(:,1)}, length(req.sample), 1) , req.sample, 'UniformOutput',0))  ); 
% subjdata= [subjdata(whichsubs,:); subjdata(end,:)];  logg.subjects=  subjdata(1:end-1,1);   logg.n_subjs= length(logg.subjects);
end
 

%% [0 Simple ad-hoc hard-coded figs]

do_simplots=0;
for o=1:1  % Simple plots
    % Calculate
    if do_simplots
        
        for o=1:1  % Archive
            do_this=0;
            if do_this
                d_plot= cell(3,3); % col: chobest, conchoPL, conchoPP. row: all trials, ncf only, cf only
                for s=1:logg.n_subjs
                    ws.d= subjdata{s,  request.session+1};
                    %         ws.d= sortrows(ws.d, [col.conchoPL col.Choice]);
                    %         ws.d(:, [col.conchoPL col.Choice col.PL1 col.PL2 col.PP1 col.PP2])
                    
                    ws.d= sortrows(ws.d, col.Trialcf);
                    %
                    % imagescnan(ws.d(:, [col.conchoPL col.conchoPP col.conchoPLPP]), 'nancolor',[0 0 0])
                    % error
                    
                    d_plot{1,1}(s)= nanmean(ws.d(:, col.conchoPLPP)==1);
                    d_plot{1,2}(s)= nanmean(ws.d(:, col.conchoPL)==1);
                    d_plot{1,3}(s)= nanmean(ws.d(:, col.conchoPP)==1);
                    
                    % No conflict
                    ws.d= d_cc{s,  request.session+1};    % Leave in draws
                    d_plot{2,1}(s)= nanmean(ws.d(:, col.conchoPLPP)==1);
                    d_plot{2,2}(s)= nanmean(ws.d(:, col.conchoPL)==1);
                    d_plot{2,3}(s)= nanmean(ws.d(:, col.conchoPP)==1);
                    
                    % Conflict trials
                    ws.d= d_ic{s,  request.session+1};
                    d_plot{3,1}(s)=   nanmean(ws.d(:, col.conchoPLPP)==1);
                    d_plot{3,2}(s)=  nanmean(ws.d(:, col.conchoPL)==1);
                    d_plot{3,3}(s)= nanmean(ws.d(:, col.conchoPP)==1);
                    
                    % cF >  no cf
                    d_plot{4,1}(s)= d_plot{3,1}(s) - d_plot{2,1}(s);
                    d_plot{4,2}(s)= d_plot{3,2}(s) - d_plot{2,2}(s);
                    d_plot{4,3}(s)= d_plot{3,3}(s) - d_plot{2,3}(s);
                end
                figure('color','w');   subplot(2,1,1); barwitherr(cellfun(@(d)std(d)./sqrt(logg.n_subjs),d_plot(1:3,:)),  cellfun(@(d)mean(d),d_plot(1:3,:)))
                %     scatter(  [1.*ones(logg.n_subjs,1); 2.*ones(logg.n_subjs,1); 3.*ones(logg.n_subjs,1)], [d_plot{4,1} d_plot{4,2} d_plot{4,3}])
                hold on, plot(0:4, [0.5 0.5 0.5 0.5 0.5], 'color',[0 0 0]), ylim([0.2 1]); ylabel('% Chose best','FontSize',20), set(gca, 'xticklabel', {'All' 'CC' 'IC'},'FontSize',20 ), legend({'PLPP' 'PL' 'PP'},'FontSize',20)
                subplot(2,1,2);  barwitherr(cellfun(@(d)std(d)./sqrt(logg.n_subjs),d_plot(4,:)),  cellfun(@(d)mean(d),d_plot(4,:)), 'y'); hold on; scatter(  [1.*ones(logg.n_subjs,1); 2.*ones(logg.n_subjs,1); 3.*ones(logg.n_subjs,1)], [d_plot{4,1} d_plot{4,2} d_plot{4,3}]);  ylabel('% Chose best','FontSize',20), set(gca, 'xticklabel', {'PLPP' 'PL' 'PP'},'FontSize',20 );  title('% on IC trial > % on CC trials'), ylim auto
                [wr.anova]=teg_repeated_measures_ANOVA([d_plot{2,2}' d_plot{2,3}' d_plot{3,2}' d_plot{3,3}'] ,[2 2], {'Conflict' 'PLPP'});   disp('ANOVA:'),  disp(wr.anova.R(:,1:4));% 	Row=Var1, Var2, TxC; Col=F, df1, df2, p
            end
            
            do_this=0;
            if do_this % No. draws?
                for s=1:logg.n_subjs
                    ws.d = d_dr{s,request.session+1};
                    dd(s,1:2) =[sum(ws.d(:,  col.PLbest)==0) sum(ws.d(:,  col.PPbest)==0)];
                end
                figure('color','w'); barwitherr(std(dd)./sqrt(logg.n_subjs), mean(dd), 'y'); set(gca,'xticklabel',{'PL', 'PP'}),  xlim([0.5 2.5]), ylabel('N draws')
                [h p ci st]=ttest(dd(:,1)-dd(:,2));%  p, st.tstat;
            end
            
            % Does PL/PP correspond more to overall PLPP?
            do_PLPP_correspond=0;
            if do_PLPP_correspond  % This need more thinking out
                d= subjdata{logg.n_subjs+1, request.session+1};  d= [d(:, [col.PL1 col.PP1 col.PLPP1]); d(:, [col.PL2 col.PP2 col.PLPP2])];
                close all,    figure('color','w'); scatter( d(:, 1),  d(:,3), 'marker', '+', 'markeredgecolor','r'), hold on, scatter( d(:, 2),  d(:,3), 'marker', 'o', 'markeredgecolor','b')
                lsline, xlabel('PL (red), PP (blue)','FontSize',20), ylabel('PLPP','FontSize',20), set(gca,'FontSize',20)
                [r p]=corr(d); disp(['PL & PLPP: r=' num2str(r(3,1)) ', p=' num2str(p(3,1)) '    '  'PP & PLPP: r=' num2str(r(3,2)) ', p=' num2str(p(3,2))])
                d= subjdata{logg.n_subjs+1, request.session+1}; [r p]=corr(d(:, [col.PLbest col.PPbest col.PLPPbest]) );
            end
            
            do_this=0;
            if do_this
                
                
                % RT slower for IC than CC?
                d_rt=nan(logg.n_subjs,4);  % CC, IC, IC-choPL, IC-choPP
                for s=1:logg.n_subjs
                    ws.d=subjdata{s, request.session+1};
                    d_rt(s,1)=  nanmean(ws.d(ws.d(:, col.Trialcf)==0, col.ChoiceRT));  % CC
                    d_rt(s,2)=  nanmean(ws.d(ws.d(:, col.Trialcf)==1, col.ChoiceRT));  % IC
                    ws.dic= ws.d(ws.d(:, col.Trialcf)==1, :);
                    d_rt(s,3)=  nanmean(ws.dic(ws.dic(:, col.PLbest)==ws.dic(:, col.Choice), col.ChoiceRT));
                    d_rt(s,4)=  nanmean(ws.dic(ws.dic(:, col.PPbest)==ws.dic(:, col.Choice), col.ChoiceRT));
                end
                figure('color','w','name', 'Mean RTs');  barwitherr(nanstd(d_rt)./sqrt(logg.n_subjs-sum(isnan(d_rt))), nanmean(d_rt), 'y'); hold on, scatter(sortrows(repmat((1:4)',logg.n_subjs, 1) ), d_rt(:))
                %         ylim([1.2 2]),
                ylabel('Mean z-scored RT','FontSize',20); title('Mean RTs by trial type','FontSize',20);  set(gca,'FontSize',20,'xticklabel',{'CC';'IC';'IC-PL';'IC-PP'})
                disp('t-tests:   ------')
                [h p ci st]=ttest(d_rt(:,1), d_rt(:,2));  disp(['CC vs IC: t(' num2str(st.df) ')= ' num2str(st.tstat)  ' , p='  num2str(p)])
                [h p ci st]=ttest(d_rt(:,1), d_rt(:,3));  disp(['CC vs IC-PL: t(' num2str(st.df) ')= ' num2str(st.tstat)  ' , p='  num2str(p)])
                [h p ci st]=ttest(d_rt(:,1), d_rt(:,4));  disp(['CC vs IC-PP: t(' num2str(st.df) ')= ' num2str(st.tstat)  ' , p='  num2str(p)])
                [h p ci st]=ttest(d_rt(:,3), d_rt(:,4));  disp(['IC-PL vs IC-PP: t(' num2str(st.df) ')= ' num2str(st.tstat)  ' , p='  num2str(p)])
                disp('NP tests:   ------')
                [p h st]=signrank(d_rt(:,1)-d_rt(:,2));   disp(['CC vs IC: z= ' num2str(st.zval)  ' , p='  num2str(p)])
                [p h st]=signrank(d_rt(:,1)-d_rt(:,3));  disp(['CC vs IC-PL: z= ' num2str(st.zval)  ' , p='  num2str(p)])
                [p h st]=signrank(d_rt(:,1)-d_rt(:,4));  disp(['CC vs IC-PP: z= ' num2str(st.zval)  ' , p='  num2str(p)])
                [p h st]=signrank(d_rt(:,3)-d_rt(:,4));  disp(['IC-PL vs IC-PP: z= ' num2str(st.zval)  ' , p='  num2str(p)])
            end
            
        end
        
        % Whats the relationship between PL and PP?
        do_plpp_rel=0;
        if do_plpp_rel
            close all hidden;   d_plpp=cell(logg.n_subjs,2);     w.bins= unique(subjdata{logg.n_subjs+1, request.session+1}(:, [col.PL1 col.PP1 col.PL2 col.PP2]));  d_plpp{logg.n_subjs+1,2}=zeros(length(w.bins));
            for s=1:logg.n_subjs
                ws.d= subjdata{s, request.session+1}; 
                ws.dd= [ws.d(:, [col.PL1 col.PP1 col.PLPP1 ]);  ws.d(:, [col.PL2 col.PP2 col.PLPP2])];
                d_plpp{s,1}= ws.dd; d_plpp{s,2}=nan(length(w.bins));
                
                for i1=1:length(w.bins)
                    for i2=1:length(w.bins)
                        d_plpp{s,2}(i1, i2)= sum( ws.dd(:, 1)==w.bins(i1) & ws.dd(:, 2)==w.bins(i2) );
                    end
                end
                
                d_plpp{logg.n_subjs+1,1}= [d_plpp{logg.n_subjs+1,1}; [ws.dd ones(size(ws.dd,1),1)*s]  ];
                d_plpp{logg.n_subjs+1,2}=d_plpp{logg.n_subjs+1,2}+d_plpp{s,2};
            end
            d_plpp{logg.n_subjs+1,2}=d_plpp{logg.n_subjs+1,2}./logg.n_subjs;
            figure('color','w'); imagescnan(d_plpp{logg.n_subjs+1,2}), colorbar, axis square , set(gca,'FontSize',20)
            figure('color','w');   for i=1:8
                subplot(2,4,i)
                imagescnan(d_plpp{logg.n_subjs+1,2}), colorbar, axis square , set(gca,'FontSize',20)
                title(['Plot re-scaled: 0-' num2str(12-i)],'FontSize',20)
                caxis([0 12-i])
            end
            
            % Straightforward correlations between PLPP
            r_plppcorr=  sortrows([logg.subjects num2cell(cellfun(@(x)corr(x(:,1),x(:,2)),  d_plpp(1:logg.n_subjs,1)))],2);
            f.plotcols= 5;  f.figwidth= 800; f.figheight=500; f.fontsize=10; f.fontsize_title=30;f.fontname='PT Sans Caption';
            f.subplot_VerHorz=[0.08 0.02]; f.fig_BotTop=[0.1 0.1]; f.fig_LeftRight=[0.15 0.15];
            figure('Name', 'PL x PP, ntrials', 'Position', [100 550 f.figwidth f.figheight], 'Color', 'w');
            for ss=1:logg.n_subjs
                s= find(strcmp(logg.subjects, r_plppcorr{ss,1}));
                subtightplot(ceil(logg.n_subjs/f.plotcols),  f.plotcols, ss ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);
                %                 imagescnan(d_plpp{s,2}, 'nancolor',[0 0 0]), colorbar, axis square, axis off
                %         caxis([min(d_subntrials{logg.n_subjs,iv}(:)) max(d_subntrials{logg.n_subjs,iv}(:))])
                scatter( d_plpp{s,1}(:,1), d_plpp{s,1}(:,2)); h=lsline; set(h,'color',[0 0 0],'LineWidth',10)
                xlim([1 10]),ylim([1 10]), axis square
                title(logg.subjects{s},'FontSize', f.fontsize)
            end
            disp(['Average r PL&PP:  mean= '  num2str(mean(cell2mat(r_plppcorr(:,2)))) ', SD=' num2str(std(cell2mat(r_plppcorr(:,2))))]); % figure('color','w');  barwitherr( std(cell2mat(r_plppcorr(:,2)))./sqrt(logg.n_subjs), mean(cell2mat(r_plppcorr(:,2))), 'y'); xlim([0.5 1.5]); set(gca,'Fontsize',20); ylabel('PL & PP r','FontSize',20)
            %            r_plppcorr(1:10,1)
            %             r_plppcorr(11:20,1)
            error('done')
        end
        
        
        % Whats the relationship between PL/PP and cf?
        do_cfraw_rel=1;
        if do_cfraw_rel
            close all hidden;   d_score=cell(logg.n_subjs,2);     
%             w.bins= unique(subjdata{logg.n_subjs+1, request.session+1}(:, [col.PL1 col.PP1 col.PL2 col.PP2]));  d_score{logg.n_subjs+1,2}=zeros(length(w.bins));
            request.scores={
%                 'PLbincf';'PLbinscore';'choPLbin';  'PPbincf'; 'PPbinscore';'choPPbin';
                'PLbincf'; 'choPLbin';  'PPbincf'; 'choPPbin';
%                 'PLcf';'PLscore';'choPL';  'PPcf'; 'PPscore';'choPP';
                };   
            col.PLscore=structmax(col)+1;  col.PPscore=structmax(col)+1;  col.PLbinscore=structmax(col)+1;  col.PPbinscore=structmax(col)+1; 
            request.scols= cell2mat(cellfun(@(x,y)eval(['y.' x]), request.scores,  repmat({col}, length(request.scores),1), 'UniformOutput',0));  d_score{logg.n_subjs+1,1}=[]; 
            for s=1:logg.n_subjs  % d_scores: raw, r matrix, p matrix
                ws.dd= subjdata{s, request.session+1}; 
                ws.d=[ws.dd; ws.dd];  ws.d(:, [col.PLscore col.PPscore])= [ws.dd(:, [col.PL1 col.PP1]); ws.dd(:, [col.PL2 col.PP2])];
                ws.d=[ws.dd; ws.dd];  ws.d(:, [col.PLbinscore col.PPbinscore])= [ws.dd(:, [col.PLbin1 col.PPbin1]); ws.dd(:, [col.PLbin2 col.PPbin2])];
                d_score{s,1}=  ws.d(:,  request.scols);
                [d_score{s,2} d_score{s,3}]= corr(d_score{s,1});
                d_score{s,2}=d_score{s,2}(:); d_score{s,2}(d_score{s,3}(:)>0.1)=nan; 
                d_score{s,2}=reshape(d_score{s,2}, size(d_score{s,3}));
                
                d_score{logg.n_subjs+1,1}= [d_score{logg.n_subjs+1,1}; [d_score{s,1} ones(size(d_score{s,1},1),1)*s]  ];
            end 
            [d_score{logg.n_subjs+1,2} d_score{logg.n_subjs+1,3}]= corr(d_score{logg.n_subjs+1,1}(:, 1:end-1)); 
            d_score{logg.n_subjs+1,2}= d_score{logg.n_subjs+1,2}(:);  d_score{logg.n_subjs+1,2}(d_score{logg.n_subjs+1,3}(:)>0.1)=nan; 
            d_score{logg.n_subjs+1,2}= reshape(d_score{logg.n_subjs+1,2}, size(d_score{logg.n_subjs+1,3}));
            figure('color','w');   subplot(1,2,1),  imagescnan(d_score{logg.n_subjs+1,2}), colorbar, axis square , set(gca,'FontSize',20);set(gca,  'xtick', 1:length(request.scores),  'xticklabel', 1:length(request.scores),   'ytick', 1:length(request.scores),  'yticklabel',  cellfun(@(x, y)[x '  [' y ']'], request.scores ,  num2cell(num2str((1:length(request.scores))')), 'UniformOutput',0), 'FontSize',15); title('Mean r','FontSize',20)
            subplot(1,2,2),  imagescnan(d_score{logg.n_subjs+1,3}), colorbar, axis square , set(gca,'FontSize',20); set(gca,  'xtick', 1:length(request.scores),  'xticklabel', 1:length(request.scores),   'ytick', 1:length(request.scores),  'yticklabel',  cellfun(@(x, y)[x '  [' y ']'], request.scores ,  num2cell(num2str((1:length(request.scores))')), 'UniformOutput',0), 'FontSize',15); title('p','FontSize',20); caxis([0 0.1]); 
            
            
            % Individual plots S traightforward correlations between PLPP
%             d_plot(:,1)= cellfun(@(x)x(:, [0+1 0+3]), d_score(:,1), 'UniformOutput',0);
%             d_plot(:,2)= cellfun(@(x)x(:, [3+1 3+3]), d_score(:,1), 'UniformOutput',0);
            f.plotcols= 5;  f.figwidth= 800; f.figheight=500; f.fontsize=10; f.fontsize_title=30;f.fontname='PT Sans Caption';
            f.subplot_VerHorz=[0.03 0.05]; f.fig_BotTop=[0.05 0.05]; f.fig_LeftRight=[0.15 0.05];
            figure('Name', 'Score correlation', 'Position', [100 550 f.figwidth f.figheight], 'Color', 'w');
            for s=1:logg.n_subjs
%                 s= find(strcmp(logg.subjects, r_corr{ss,1}));
                subtightplot(ceil(logg.n_subjs/f.plotcols),  f.plotcols, s ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);
                %                 imagescnan(d_score{s,2}, 'nancolor',[0 0 0]), colorbar, axis square, axis off
                %         caxis([min(d_subntrials{logg.n_subjs,iv}(:)) max(d_subntrials{logg.n_subjs,iv}(:))])
                 imagescnan(d_score{s,2}), colorbar, axis square
                 set(gca,  'xtick', 1:length(request.scores),  'ytick', 1:length(request.scores)) 
                 if rem(s,f.plotcols)==1;  set(gca,  'yticklabel',  cellfun(@(x, y)[x '  [' y ']'], request.scores ,  num2cell(num2str((1:length(request.scores))')), 'UniformOutput',0), 'FontSize',15);  end
                 title(logg.subjects{s},'FontSize', f.fontsize)
            end
            %            r_corr(1:10,1)
            %             r_corr(11:20,1)
            error('done')
        end
        
        
        
    end
end

%% [1a: unidimensional IV] Range of variables experienced across task (mostly IVs)

for o1=1:1
    dop_vec_vars=0;
    
    if dop_vec_vars
        request.var ={
            'PLcf' 'PL cf'; 'PPcf' 'PP cf';  % 'PLPPcf' 'PLPP cf';                   % Raw scores
            %         'PLPPcf' 'PLPP cf';
            %         'ad1m2.PL' 'absdiff PLbin'; 'ad1m2.PP' 'absdiff PPbin'     % Binned scores
            %         'PLbincf' 'PLbin cf'; 'PPbincf' 'PPbin cf'
            %         'PLbinvcf' 'PLbin vcf'; 'PPbinvcf' 'PPbin vcf'
            %         'PLvcf' 'PL vcf'; 'PPvcf' 'PP vcf';                  % Conflict scaled by valu e
            %         'PLPPvcf' 'PLPP vcf';
            %         'ad1m2.PL' 'absdiff PL'; 'ad1m2.PP' 'absdiff PP';   % Absolute difference
            %         'ad1m2.PLPP' 'absdiff PLPP';
            %         'choPL' 'cho PL'; 'unchoPL' 'uncho PL'; 'choPP' 'cho PP'; 'unchoPP' 'uncho PP';  % Cho/uncho
            %         'choPLPP' 'cho PLPP'; 'unchoPLPP' 'uncho PLPP';
            
            };
        request.nvar=size(request.var,1);  request.varcols=  cell2mat(cellfun(@(x,col)eval(['col.'  x]),  request.var(:,1) ,  repmat({col}, size(request.var,1),1), 'UniformOutput',0) );
        d_var(logg.n_subjs+1,1:3)= cell(1,3) ;   % Variables to look at. Cols = order in request.var
        for s=1:logg.n_subjs
            ws.d= subjdata{s,request.session+1};
            d_var{s,1}= ws.d(:, request.varcols);
            d_var{s,2}= mean(d_var{s,1});
            d_var{s,3}= std(d_var{s,1});
            d_var{logg.n_subjs+1,1}=[d_var{logg.n_subjs+1,1}; d_var{s,1}];
            d_var{logg.n_subjs+1,2}=[d_var{logg.n_subjs+1,2}; d_var{s,2}];
            d_var{logg.n_subjs+1,3}=[d_var{logg.n_subjs+1,3}; d_var{s,3}];
            ws=[];
        end
        % d_var{logg.n_subjs+1,1}= mean(d_var{logg.n_subjs+1,1});
        % d_var{logg.n_subjs+1,2}= mean(d_var{logg.n_subjs+1,2});
        % d_var{logg.n_subjs+1,3}= mean(d_var{logg.n_subjs+1,3});
        
        figure('color','w')
        barwitherr(std(d_var{logg.n_subjs+1,2})./sqrt(logg.n_subjs),  mean(d_var{logg.n_subjs+1,2}))
        set(gca, 'xticklabel', request.var(:,2)), title('Mean IV values experienced by subjects', 'FontSize', 20)
    end
    
    
    %% [1b: unidimensional IV] Plot DVs against PL & PP scores
    
    % isthere=zeros(10,20); for s=1:20; %     isthere(unique(subjdata{s,3}(:, col.PPcf)),s)=1; % end; sum(isthere,2)
    
    dop_vec_ivdv=0;
    if dop_vec_ivdv
        col.pChoPLPP_PL = col.pChoPLPP; col.pChoPLPP_PP = col.pChoPLPP;
        col.nTrialsPL= col.nTrials; col.nTrialsPP= col.nTrials;
        request.ivs= {
            %         'mar'; 'marbin';
            'cf';
            %         'bincf'
            };     % name, cols, iv bins (PL, PP). Create binned IVs in preproc if you want to plot binned
        request.nbins=[7 10 3]; % Empty to auto
        request.dvs= {
            %          'pChoPLPP_';  %  'pCho';
            'RTcho'
            }; % pCho is choice-variant.
        
        d_plotsubd=cell(size(request.ivs,1), 2, 2, logg.n_subjs);   % {iv, iv_PL/PP, bin, subject}
        d_plotd= nan(size(request.ivs,1), 2, 2, size(request.dvs,1)*2, logg.n_subjs);   % (iv, iv_PL/PP, bin, dv, s);
        
        for o=1:1  % Setup
            % Columns
            request.dv_cols =   cell2mat( cellfun(@(x, col)[eval(['col.' x 'PL']) eval(['col.' x 'PP'])], request.dvs,  repmat({col},size(request.dvs,1),1), 'UniformOutput',0));
            request.iv_cols =  cell2mat( cellfun(@(x, col)[eval(['col.PL' x]) eval(['col.PP' x])], request.ivs ,  repmat({col},  size(request.ivs,1),1),  'UniformOutput',0) );
            %     request.iv_cols =  cellfun(@(x, col)[eval(['col.PL' x]) eval(['col.PP' x])], request.ivs(strcmp(request.ivs(:,1), 'ad1m2')==0),  repmat({col},  sum(strcmp(request.ivs(:,1), 'ad1m2')==0),1), 'UniformOutput',0);  request.ivs(find(strcmp(request.ivs, 'ad1m2')==0),2) = request.iv_cols;
            %     request.iv_cols =  cellfun(@(x, col)[eval(['col.' x '.PLbin']) eval(['col.' x '.PPbin'])], request.ivs(strcmp(request.ivs(:,1), 'ad1m2')==1), repmat({col},  sum(strcmp(request.ivs(:,1), 'ad1m2')==1),1), 'UniformOutput',0);  request.ivs(find(strcmp(request.ivs, 'ad1m2')==1),2) = request.iv_cols;
            %     request.iv_cols=cell2mat(request.ivs(:,2));  % PL, then PP
            
            % IV bins
            for iv=1:size(request.ivs,1)
                request.ivs{iv,3}= {unique( subjdata{logg.n_subjs+1,request.session+1}(:, request.iv_cols(iv,1)) ) unique( subjdata{logg.n_subjs+1,request.session+1}(:, request.iv_cols(iv,2)) )};
            end
        end
        for s=1:logg.n_subjs
            ws.d= subjdata{s,request.session+1};
            
            % Split data by IV bins
            for iv=1:size(request.ivs,1)
                for p=1:2
                    wi.bins =unique(request.ivs{iv,3}{p});
                    for b=1:length(wi.bins )
                        d_plotsubd{iv, p,b,s}= ws.d( ws.d(:,  request.iv_cols(iv,p))==wi.bins(b), :);
                        d_plotsubd{iv, p,b,s}(:, col.nTrials)=0;
                        if isempty(d_plotsubd{iv, p,b,s})==0
                            % Calculate new values
                            wi.ntrials= size(d_plotsubd{iv, p,b,s},1);
                            wi.choPL= find(d_plotsubd{iv, p,b,s}(:, col.choPL)>d_plotsubd{iv, p,b,s}(:, col.unchoPL));  % Choice congruent w PL
                            wi.choPP= find(d_plotsubd{iv, p,b,s}(:, col.choPP)>d_plotsubd{iv, p,b,s}(:, col.unchoPP));   % Choice congruent w PP
                            wi.choPLPP= find(d_plotsubd{iv, p,b,s}(:, col.choPLPP)>d_plotsubd{iv, p,b,s}(:, col.unchoPLPP));
                            wi.choeven= find( (d_plotsubd{iv, p,b,s}(:, col.choPL)==d_plotsubd{iv, p,b,s}(:, col.unchoPL)) + (d_plotsubd{iv, p,b,s}(:, col.choPP)==d_plotsubd{iv, p,b,s}(:, col.unchoPP)));  % PL and PP are both even between cho and uncho
                            wi.choerr= find(1- ((d_plotsubd{iv, p,b,s}(:, col.choPL)>d_plotsubd{iv, p,b,s}(:, col.unchoPL)) + (d_plotsubd{iv, p,b,s}(:, col.choPP)>d_plotsubd{iv, p,b,s}(:, col.unchoPP)) + (d_plotsubd{iv, p,b,s}(:, col.choPL)==d_plotsubd{iv, p,b,s}(:, col.unchoPL)) + (d_plotsubd{iv, p,b,s}(:, col.choPP)==d_plotsubd{iv, p,b,s}(:, col.unchoPP))) );  % Apparent errors: uncho better than cho on both dimensions
                            %                         if length([wi.choPL ;wi.choPP; wi.choeven; wi.choerr])> wi.ntrials; disp('FLAG!'); end
                            d_plotsubd{iv, p,b,s}(:, [col.pChoPL col.pChoPP col.pChoPLPP])= repmat([length(wi.choPL)/wi.ntrials length(wi.choPP)/wi.ntrials length(wi.choPLPP)/wi.ntrials], wi.ntrials ,1);
                            d_plotsubd{iv, p,b,s}(:, [col.RTchoPL col.RTchoPP])=nan;
                            
                            d_plotsubd{iv, p,b,s}(:, col.nTrials)=wi.ntrials;
                            
                            d_plotsubd{iv, p,b,s}(wi.choPL, col.RTchoPL)=  d_plotsubd{iv, p,b,s}(wi.choPL, col.ChoiceRT);
                            d_plotsubd{iv, p,b,s}(wi.choPP, col.RTchoPP)=  d_plotsubd{iv, p,b,s}(wi.choPP, col.ChoiceRT);
                            ws.dmat{iv, p,b,s}= nanmean(d_plotsubd{iv, p,b,s},1);
                            
                            for dv=1:size(request.dvs,1)  % Auto-readout DVs
                                d_plotd(iv, p,b, dv,s) =  ws.dmat{iv, p,b,s}(:, request.dv_cols(dv,p));
                            end
                        end
                        
                    end
                    wi=[];
                end
            end
            
            ws=[];
        end
        d_plotm = nanmean(d_plotd,5);    % nanmean: ignore how many subjects have good data
        
        % plot
        f.plotcols=length(request.dvs);  f.figwidth= 500; f.figheight=1400; f.fontsize=20; f.fontsize_title=30;f.fontname='PT Sans Caption';
        f.subplot_VerHorz=[0.15 0.15]; f.fig_BotTop=[0.15 0.15]; f.fig_LeftRight=[0.1 0.1];
        figure('Name', ' ', 'Position', [100 550 f.figwidth f.figheight], 'Color', 'w');  k=1;
        for iv=1:size(request.ivs,1)
            for dv=1:size(request.dvs,1)
                subtightplot(size(request.ivs,1),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
                imagescnan([squeeze(d_plotm(iv, 1, :,dv)) squeeze(d_plotm(iv, 2, :,dv))], 'nancolor',[0 0.6 0]), colorbar
                colormap 'pink'
                set(gca, 'xticklabel', {'PL' 'PP'},'Fontsize',f.fontsize), title([request.ivs{iv,1} ' - ' request.dvs{dv,1}])
                
            end
        end
    end
    
end

%% [2] Plot requested x for PLdimension vs PP dimension
% Always PL first, PP second 

col.PLcho=col.choPL;  col.PPcho=col.choPP;  col.PLuncho=col.unchoPL;  col.PPuncho=col.unchoPP;  col.PLPPcho=col.choPLPP;    col.PLPPuncho=col.unchoPLPP;
col.PLchobin=col.choPLbin;  col.PPchobin=col.choPPbin;  col.PLunchobin=col.unchoPLbin;  col.PPunchobin=col.unchoPPbin;  col.PLPPchobin=col.choPLPPbin;    col.PLPPunchobin=col.unchoPLPPbin;
col.cfweightPL= col.cfweight;  col.cfweightPP= col.cfweight;  col.cfweightbinPL= col.cfweightbin;  col.cfweightbinPP= col.cfweightbin; 
col.colmaxf=structmax(col);    
col.nTrialschoPL= col.colmaxf+1;  col.nTrialschoPP= col.colmaxf+2;  col.pTrialschoPL= col.colmaxf+3;  col.pTrialschoPP= col.colmaxf+4; 

close all hidden
dop_plpp_ivdv=1;
if dop_plpp_ivdv
    request.ivs={
        'bincf'; 
%         'cf';
%             'mar';  
% 'cho';
% 'uncho'    
%         'marbin'; 
        'chobin'; 
% 'unchobin'
%         'binvcf'; 'vcf'
        };
    
%     request.ivs={'bincf';   'cho'; 'chobin';  'marbin'; }; 
    request.dvs =  {  % DVs: If any of these need to be calculated on the fly (i.e. within specified IV cells), make sure they are computed below!
%             'pChoPLPP';  'pChoPL'; 'pChoPP';
        'pChoPLPPbin'; 'pChoPLbin'; 'pChoPPbin';
        'ChoiceRT';
% 'RTchoPL'; 'RTchoPP'
%         'nTrials'; 'nTrialschoPL'; 'nTrialschoPP';
%         'pTrialschoPL'; 'pTrialschoPP';
%         'PLcf';'PPcf';
        } ;
    
    d_subd=repmat({cell(2,2)}, logg.n_subjs, size(request.ivs,1));   % {subject, iv}{PLbin, PPbin}: all data for that cell (for subject)
    d_subntrials=repmat({nan(2,2)}, logg.n_subjs+1, size(request.ivs,1));   % {subject, iv}{PLbin, PPbin}: n trials for that cell (for subject)
    d_dv=repmat({cell(2,2, logg.n_subjs )}, size(request.ivs,1),   size(request.dvs,1));   % {iv, dv}{PLbin, PPbin, subject}: all data for that cell (for subject)
    for o=1:1  % Setup
        % [Columns] IVs are assumed to have PL & PP dimensions; DVs are assumed to not.
        request.iv_cols =  cell2mat( cellfun(@(x, col)[eval(['col.PL' x]) eval(['col.PP' x])], request.ivs ,  repmat({col},  size(request.ivs,1),1),  'UniformOutput',0) );
        request.dv_cols =   cell2mat( cellfun(@(x, col)eval(['col.' x]), request.dvs,  repmat({col},size(request.dvs,1),1), 'UniformOutput',0));
        
        % IV bins
        for iv=1:size(request.ivs,1)
            request.ivs{iv,3}= {unique( subjdata{logg.n_subjs+1,request.session+1}(:, request.iv_cols(iv,1)) ) unique( subjdata{logg.n_subjs+1,request.session+1}(:, request.iv_cols(iv,2)) )};
            request.ivs{iv,3}{1} = request.ivs{iv,3}{1}(find(1-isnan(request.ivs{iv,3}{1})));
            request.ivs{iv,3}{2}=request.ivs{iv,3}{2}(find(1-isnan(request.ivs{iv,3}{2})));
            d_subntrials{logg.n_subjs+1, iv}= zeros(length(request.ivs{iv,3}{1}), length(request.ivs{iv,3}{2}));
        end
        
    end
    for s=1:logg.n_subjs % Extract data 
        ws.d= subjdata{s,request.session+1};
        %     ws.d=sortrows(ws.d, [col.PPcf col.PLcf]);

        
        % Bin raw data into IV cells + calculate new quantities if you wan
        for iv=1:size(request.ivs,1)
            for ipl=1:length(request.ivs{iv,3}{1})
                for ipp=1:length(request.ivs{iv,3}{2})
                    wc.d= ws.d(   ws.d(:,  request.iv_cols(iv, 1)) == request.ivs{iv,3}{1}(ipl) & ws.d(:,  request.iv_cols(iv, 2)) == request.ivs{iv,3}{2}(ipp) , :);
                    d_subntrials{s,iv}(ipl, ipp)=size(wc.d,1);
                    d_subntrials{logg.n_subjs+1,iv}(ipl, ipp)=d_subntrials{logg.n_subjs+1,iv}(ipl, ipp)+  d_subntrials{s,iv}(ipl, ipp);
                    wc.d(:, col.nTrials)=size(wc.d,1);
                    
%                     
%         if s==5 && ipl==10 & ipl==10
%             
%             ; error; 
%         
%                          wc.d([col.PL1 col.PP1 col.PL2 col.PP2 col.Choice  col.choPL col.unchoPL col.choPP col.unchoPP])
%                         
%         end 
                    
                    % Hard-code-compute new quantities, within cell
                    %  Quantities that are coded within trial should be completed earlier
                    if isempty(wc.d)==0
                        wc.d(:, col.pChoPL)= sum(wc.d(:, col.choPL) > wc.d(:, col.unchoPL))/size(wc.d,1);
                        wc.d(:, col.pChoPP)= sum(wc.d(:, col.choPP) > wc.d(:, col.unchoPP))/size(wc.d,1);
                        wc.d(:, col.pChoPLPP)= sum(wc.d(:, col.choPLPP) > wc.d(:, col.unchoPLPP))/size(wc.d,1);

                        
                        
                        % if size(wc.d,1)>5; error; end
                        wc.d(:, col.pChoPLbin)= sum(wc.d(:, col.choPLbin) > wc.d(:, col.unchoPLbin))/size(wc.d,1);
                        wc.d(:, col.pChoPPbin)= sum(wc.d(:, col.choPPbin) > wc.d(:, col.unchoPPbin))/size(wc.d,1);
                        wc.d(:, col.pChoPLPPbin)= sum(wc.d(:, col.choPLPPbin) > wc.d(:, col.unchoPLPPbin))/size(wc.d,1);
                        
                        
%                         wc.d(wc.d(:, col.Choice)== wc.d(:, col.PLbest), col.nTrialschoPL  )=sum(wc.d(:, col.Choice)== wc.d(:, col.PLbest)); 
%                         wc.d(wc.d(:, col.Choice)== wc.d(:, col.PPbest), col.nTrialschoPP  )=sum(wc.d(:, col.Choice)== wc.d(:, col.PPbest)); 
%                         wc.d(wc.d(:, col.Choice)== wc.d(:, col.PLbest), col.pTrialschoPL  )=mean(wc.d(:, col.Choice)== wc.d(:, col.PLbest)); 
%                         wc.d(wc.d(:, col.Choice)== wc.d(:, col.PPbest), col.pTrialschoPP  )=mean(wc.d(:, col.Choice)== wc.d(:, col.PPbest)); 
                        
                        
                        
                        d_subd{s,iv}{ipl, ipp} = wc.d;
                    else d_subd{s,iv}{ipl, ipp} = nan;  % If you change this, it will add zeros to all unfilled columns 
                    end
                    
                end
            end
        end
        
        % Auto-read DVs
        for iv=1:size(request.ivs,1)
            for ipl=1:length(request.ivs{iv,3}{1})
                for ipp=1:length(request.ivs{iv,3}{2})
                    for dv=1:size(request.dvs,1)  % Auto-read DVs for non-empty cells (no on-the-fly calculation
                        if isempty(d_subd{s,iv}{ipl, ipp} )==1 | (isnan(d_subd{s,iv}{ipl, ipp}) & size(d_subd{s,iv}{ipl, ipp},2)==1)   % empty d_subd cells had a nan inserted
                            d_dv{iv,dv}{ipl, ipp, s} =nan;
                        else d_dv{iv,dv}{ipl, ipp, s} = nanmean(d_subd{s,iv}{ipl, ipp}(:, request.dv_cols(dv)));
                        end
                    end
                end
            end
        end
        
        ws=[];
    end
    
%     close all hidden
    % Plot
%     f.axis=[1.5 1.9];   % Comment out to omit 
    f.plotcols=  size(request.dvs,1); % f.plotcols= 4; 
    f.figwidth= 800; f.figheight=500; f.fontsize=13; f.fontsize_title=30;f.fontname='PT Sans Caption';
%     f.subplot_VerHorz=[0.15 0.05]; f.fig_BotTop=[0.15 0.1]; f.fig_LeftRight=[0.2 0.2];    
    f.subplot_VerHorz=[0.15 0.01]; f.fig_BotTop=[0.05 0.05]; f.fig_LeftRight=[0.25 0.25];
    figure('Name', 'PLxPP', 'Position', [100 550 f.figwidth f.figheight], 'Color', 'w'); k=1;
    for iv=1:size(request.ivs,1)
        for dv=1:size(request.dvs,1)
            wc.d = d_dv{iv,dv};
            wc.plot = nanmean(cell2mat(wc.d) ,3); % This is what you want 
%             disp([request.ivs{iv,1} '   '  request.dvs{dv,1}]), disp(wc.plot)
%             wc.plot=wc.plot (2:end, 2:end); disp('Hard code ranges') 
            
            if isempty(strfind(request.dvs{dv}, 'RT'))==0; wc.nancol=[0 0.4 0]; else wc.nancol=[0 0.0 0]; end
            subtightplot(size(request.ivs,1), f.plotcols , k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
            imagescnan(wc.plot  , 'nancolor', wc.nancol), colorbar, axis square
            title(['[' request.ivs{iv} ']   ' request.dvs{dv}],'FontSize', f.fontsize)
            set(gca,'FontSize', f.fontsize)
            %         ylabel(['PL ' request.ivs{iv}],'FontSize', f.fontsize), xlabel(['PP ' request.ivs{iv}],'FontSize', f.fontsize)
            
            set(gca, 'ytick',  1:length(request.ivs{iv,3}{1}), 'yticklabel', num2cell(request.ivs{iv,3}{1}),  'xtick',  1:length(request.ivs{iv,3}{2}),  'xticklabel', num2cell(request.ivs{iv,3}{2}))
            
%             if isempty(strfind(request.dvs{dv}, 'RT'))==0; colormap 'pink'; end
            if isfield(f,'axis'),  caxis(f.axis),    disp('Artificial cxis!'),end
            % eror
        end
    end
    
    
%     close all hidden
    error('done')
    
    
% % [ MISCELLANY FROM THIS SECTION ] ##########
%     close all
    % Plot individual
    do_plotindiv= 1;
    if do_plotindiv,   iv=1; dv=2 ;f.plotcols= 5; f.figwidth= 800; f.figheight=250; f.fontsize=15; f.fontsize_title=20;f.fontname='PT Sans Caption';
        disp('wide format'); f.subplot_VerHorz=[0.05 0.01]; f.fig_BotTop=[0.05 0.05]; f.fig_LeftRight=[0.15 0.25];
%         disp('superwide format'); f.subplot_VerHorz=[0.01 0.01]; f.fig_BotTop=[0.01 0.01]; f.fig_LeftRight=[0.05 0.05]; f.plotcols=10; f.figwidth= 1000; f.figheight=200;
%         disp('long format'); f.subplot_VerHorz=[0.01 0.025]; f.fig_BotTop=[0.05 0.05]; f.fig_LeftRight=[0.15 0.25]; f.plotcols= 2; f.figwidth= 250; f.figheight=800;
        figure('Name', 'PLxPP sub', 'Position', [800 0 f.figwidth f.figheight], 'Color', 'w'); k=1;
        for s=1:logg.n_subjs
            subtightplot(ceil(logg.n_subjs/f.plotcols) , f.plotcols , k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
            imagescnan(cell2mat(d_dv{iv,dv}(:,:,s)) , 'nancolor', [0 0 0]), colorbar, axis square
            %     title(logg.subjects{s},'FontSize', f.fontsize)
            axis off, set(gca,'FontSize', f.fontsize)
            %     ylabel(['PL ' request.ivs{iv}],'FontSize', f.fontsize), xlabel(['PP ' request.ivs{iv}],'FontSize', f.fontsize)
            
            % Mark which model (spacing for 3-level IVs)
            if s==1;   text(-0.25,15.4, [ 'IV: ' request.ivs{iv,1} '   DV: ' request.dvs{dv}], 'Fontsize',f.fontsize_title,'rotation',90);   end
%             if s==1;   text(-0.05,-0.8, [ 'IV: ' request.ivs{iv,1} '   DV: ' request.dvs{dv}], 'Fontsize',f.fontsize_title);   end
        end
    end
     
    % Statistics
    dostats=0;
    if dostats, disp('################### Plot statistics ###################'); d_anova=   cell(size(request.ivs,1),size(request.dvs,1)); r_anova= [[{' '}; request.ivs(:,1)] [request.dvs(:,1)';  cell(size(request.ivs,1),size(request.dvs,1))]];
        for iv=1:size(request.ivs,1)
            for dv=1:size(request.dvs,1)
                wc.d = d_dv{iv,dv};
                wc.da=[];
                for pl=1:length(request.ivs{iv,3}{1})
                    for pp=1:length(request.ivs{iv,3}{2})
                        wc.da=[wc.da  cell2mat( squeeze(wc.d(pl,pp,:)))];
                    end
                end
                d_anova{iv,dv}=wc.da;
                
                % ANOVA
                %         disp( ['PL' request.ivs{iv,1} ' x '   'PP' request.ivs{iv,1}  ' (' num2str(length(request.ivs{iv,3}{1})) 'x' num2str(length(request.ivs{iv,3}{2})) ') ANOVA on   ' request.dvs{dv,1}])
                disp( ['PL' request.ivs{iv,1} ' x PP' request.ivs{iv,1} ' fx:      '  request.dvs{dv,1} ])
                
                [wr.anova]=teg_repeated_measures_ANOVA(d_anova{iv,dv},[length(request.ivs{iv,3}{1}) length(request.ivs{iv,3}{2})], {['PL ' request.ivs{iv,1} ] ['PP ' request.ivs{iv,1} ]});  % Row=Var1, Var2, TxC; Col=F, df1, df2, p
                r_anova{iv+1, dv+1}=wr.anova.R(:,1:4);
                disp(['ME PL' request.ivs{iv,1} ':  p= ' num2str(wr.anova.R(1,4))])
                disp(['ME PP' request.ivs{iv,1} ':  p= ' num2str(wr.anova.R(2,4))])
                disp(['PLxPP ' request.ivs{iv,1} ':  p= ' num2str(wr.anova.R(3,4))])
                disp('--------------')
                
                disp(d_anova{iv,dv}'), disp('--------------')
                figure, plot( 1:length(d_anova{iv,dv}'), d_anova{iv,dv}')
                
                
            end
        end
    end
    
    dop_cellntrials=0;   % N trials plots
    if dop_cellntrials, do_indiv=0;     f.plotcols=2; f.fontsize=30; k=1; figure('color','w', 'Name','N trials');
        for iv=1:size(request.ivs,1)
            subplot( ceil(size(request.ivs,1)/ f.plotcols),f.plotcols, k)
            d_subntrials{logg.n_subjs+1,iv}=d_subntrials{logg.n_subjs+1,iv}./logg.n_subjs;
            imagesc(d_subntrials{logg.n_subjs+1,iv}), colorbar, axis square
%             xlabel(['PP ' request.ivs{iv,1}],'FontSize', f.fontsize), ylabel(['PL ' request.ivs{iv,1}],'FontSize', f.fontsize)
            set(gca,'FontSize', f.fontsize)
%             title(['[' request.ivs{iv,1} ']  Average nTrials'],'FontSize', f.fontsize)
            set(gca, 'ytick',  1:length(request.ivs{iv,3}{1}), 'yticklabel', num2cell(request.ivs{iv,3}{1}),  'xtick',  1:length(request.ivs{iv,3}{2}),  'xticklabel', num2cell(request.ivs{iv,3}{2}))            
            k=k+1;
        end
        if do_indiv
            f.plotcols= 5;  f.figwidth= 800; f.figheight=500; f.fontsize=10; f.fontsize_title=30;f.fontname='PT Sans Caption';
            f.subplot_VerHorz=[0.05 0.02]; f.fig_BotTop=[0.1 0.1]; f.fig_LeftRight=[0.05 0.05];
            for iv=1:size(request.ivs,1)  % Individual subject plots
                figure('Name', ['Ntrials - ' request.ivs{iv,1}], 'Position', [100 550 f.figwidth f.figheight], 'Color', 'w');
                for s=1:logg.n_subjs
                    subtightplot(ceil(logg.n_subjs/f.plotcols),  f.plotcols, s ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);
                    imagescnan(d_subntrials{s,iv}), colorbar, axis square, axis off
                    caxis([min(d_subntrials{logg.n_subjs,iv}(:)) max(d_subntrials{logg.n_subjs,iv}(:))])
                    title(logg.subjects{s},'FontSize', f.fontsize)
                end
            end
        end
    end
end




  
%% End


disp('########################'), disp('Preproc decisions:');  disp(' '), disp([char(logg.define(:,1)) char(repmat({'  '}, size(logg.define,1), 1)) char(logg.define(:,2))] )

% max(d_ivsfx)
% min(d_ivsfx)

 a(:,[col.PLPP(1) col.PLPP(2) col.PLPPcf])
