% GLM on choice & choice RTs:  Logistic regression for choice is p(Option 1), and difference scores are always Option1>Option2
clear all;close all hidden; clc
where.fmri_scripts ='/Users/EleanorL/Dropbox/SCRIPPS/3b PLPP fmri';  addpath(where.fmri_scripts)
for o=1:1  % Subject subsets
    logg.not_p19={'p01_YH';'p02_MI';'p03_AY';'p04_AW';'p05_CA';'p06_BT';'p07_HC';'p08_KC';'p09_KJ';'p10_YC';'p11_BI';'p12_AL';'p13_MS';'p14_SK';'p15_MM';'p16_SH';'p17_BB';'p18_WB';'p20_LZ';};
    logg.pPLunder90={'p01_YH';'p02_MI';'p03_AY';'p05_CA';'p06_BT';'p07_HC';'p09_KJ';'p10_YC';'p11_BI';'p12_AL';'p13_MS';'p14_SK';'p16_SH';'p17_BB';'p18_WB';'p20_LZ';};   % Exclude people who choose higher PL scores >90% of the time 

end

% Request specific
logg.specificsubjects={};  
% logg.specificsubjects = logg.pPLunder90; 
% logg.specificsubjects = logg.not_p19; 

logg.specificsubjects= f_subsample('beh1','PLvPP_abs_r04'); 

logg.specificsubjects= {'p01_YH';'p02_MI';'p04_AW';'p08_KC';'p09_KJ';'p11_BI';'p12_AL';'p14_SK';'p17_BB';'p18_WB';'p19_HB';};

% request.prebinned =' b3'; 
request.prebinned =[]; 
request.dataset=['All data' request.prebinned ' (31-May-2016)'];   request.choicesession=2;    % 2= fMRI, 3= in lab 
% request.dataset=['All pilot data' request.prebinned ' (04-Nov-2015)'];  request.choicesession=1;  

for o1=1:1 % General settings and specifications
request.IVs_rescale =0;  % dont   
    % Load subjects
    where.where='/Users/EleanorL/Dropbox/SCRIPPS/3 Pleasure purpose';  
    where.beh=[where.where filesep '3 Behaviour']; 
    addpath(where.where)
    if ischar(logg.specificsubjects)==1 && strcmp(logg.specificsubjects, 'AllDataOK'), logg.specificsubjects= {'p01_YH';'p02_MI';'p03_AY';'p04_AW';'p05_CA';'p06_BT';'p07_HC';'p08_KC';'p09_KJ';'p10_YC';'p11_BI';'p13_MS';'p14_SK';'p15_MM';'p16_SH';'p17_BB';'p18_WB';'p19_HB';'p20_LZ';}; end ;
    if isempty(strfind(request.dataset, 'pilot data'))==1
        [n t r]=xlsread([where.beh filesep 'datalog_plpr.xlsx']);
        [logg.subjects logg.n_subjs logg.datalog] = f_selectsubjects(r, logg.specificsubjects,r, 'All');
    else
        logg.subjects=[cellfun(@(x)['p0' num2str(x)], num2cell(1:9),'UniformOutput',0)'; cellfun(@(x)['p' num2str(x)], num2cell(10:20),'UniformOutput',0)';];
        logg.n_subjs = length(logg.subjects);
    end
    w.d=  load(['2 Data' fs request.dataset '.mat']); subjdata =  w.d.subjdata; col = w.d.col;
    [d d subjdata] = f_selectsubjects([{'Subject'}  cell(1, size(subjdata,2)-1 )  ; subjdata ], logg.subjects,[[{'Subject'}  cell(1, size(subjdata,2)-1 ) {'All'} ]; [subjdata num2cell(ones(size(subjdata,1),1))]], 'All');  subjdata=subjdata(2:end, :);
   
    
%     % Formluation 
        if isempty( w.d.log.scorebins)==1,w.cfmax= 2*10-1; 
    else w.cfmax= 2*size(log.scorebins,1)-1;
        end         
    fcf = w.d.log.define{strcmp(w.d.log.define(:,1), 'Conflict'), 3}; 
     
    
    % Misc
    request.alliv_names= {'PL1' 'PL1'; 'PL2' 'PL2';'PP1' 'PP1'; 'PP2' 'PP2';  'd1m2.PL' 'PL diff'; 'd1m2.PP' 'PP diff';    'PLcf' 'PL cf';  'PPcf' 'PP cf';   'PLPP1' 'PLPP1';   'PLPP2' 'PLPP2';  'd1m2.PLPP' 'PLPP diff';       'PLPPcf' 'PLPP cf'; 
        'er.d1m2.social'  'Social diff';  'er.ad1m2.social'  'Social absdiff'; 
        'er.d1m2.happiness' 'Happi diff';  'er.ad1m2.happiness' 'Happi absdiff';  
        'er.d1m2.vivid'  'Vivid diff';  'er.ad1m2.vivid'  'Vivid absdiff'; 
        'er.d1m2.valence'  'Valence diff'; 'er.ad1m2.valence'  'Valence absdiff'; 
        'er.d1m2.experience'  'Exper diff'; 'er.ad1m2.experience'  'Exper absdiff'; 
        'er.d1m2.arousal'  'Arousal diff';  'er.ad1m2.arousal'  'Arousal absdiff'; 
        };       
    logg.define={}; 
    
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
    disp(['No. of subjects: ' num2str(logg.n_subjs)])
    if isempty(logg.specificsubjects)==0;  if logg.n_subjs == 19 & sum(strcmp(logg.subjects, {'p01_YH';'p02_MI';'p03_AY';'p04_AW';'p05_CA';'p06_BT';'p07_HC';'p08_KC';'p09_KJ';'p10_YC';'p11_BI';'p13_MS';'p14_SK';'p15_MM';'p16_SH';'p17_BB';'p18_WB';'p19_HB';'p20_LZ';}) ) ==19; disp('   Subjects w complete data only (n=19)'); else disp('   Subset of subjects only:'); disp(logg.specificsubjects); end;   end; disp(' ')
    disp(' '); input('Hit Enter to start      ')
    disp('=======================================================')
    
end

%%

for o=1:1 % Define quantities/preprocessing
    
    f_rescale=@(rmin, rmax,x)    rmin + (rmax-rmin).*(x- min(x(:))) ./ max(x- min(x(:))); 
    
    % Note: conflict is already computed by this point. If you want to play w other formulations, 
    % give them other names. Note also that if you want to set up other IVs that are 
    % not choice invariant (i.e. dependent on chosen option's value, not just available options' values),
    % the RT analysis needs to be flipped.
  
    % [RT preprocessing]
     logg.define(size(logg.define,1)+1, :)={'RT preproc' 'Outliners (>2SD +/- mean) are nan-ed' [] };
%     frt=@(x) x-nanmean(x);  logg.define(size(logg.define,1)+1, :)={'RT preproc' 'Mean centred within subject' frt};
    frt=@(x) log(x); logg.define(size(logg.define,1)+1, :)={'RT preproc' 'Log-transformed within subject' frt};
%     frt=@(x) zscore(x);   logg.define(size(logg.define,1)+1, :)={'RT preproc' 'Zscored within subject' frt};
%     frt=@(x) 1./x;   logg.define(size(logg.define,1)+1, :)={'RT preproc' 'Reciproc' frt};
%     frt=@(x)x;   logg.define(size(logg.define,1)+1, :)={'RT preproc' 'No transform' frt};
end

% Calculate new scores  
eval(['col=col.s' num2str(request.choicesession) ';'])  % Which session? 
for o=1:1 % New scores columns 
    
    % New scores 
    col.PL1=col.PL(1); col.PL2=col.PL(2);  col.PP1=col.PP(1); col.PP2=col.PP(2); 
    col.colmax=structmax(col);   
    col.Subject=col.colmax+1;
    col.PLcf=col.colmax+2;        % Indices of conflict 
    col.PPcf=col.colmax+3;
    col.PLPPcf=col.colmax+4;
    col.PLPP1 = col.colmax+5;
    col.PLPP2 = col.colmax+6;     
    col.PLdraw=col.colmax+7; 
    col.PPdraw=col.colmax+8; 
    col.PLPPdraw=col.colmax+9; 
end
subjdata{logg.n_subjs+1,request.choicesession+1}=[];   d_cf=subjdata;d_ncf=d_cf; d_ncf_nd=d_cf;  d_ncf_drawPL=d_ncf; d_ncf_drawPP=d_ncf; d_nd=d_cf; d_dr=d_cf; d_rt=cell(logg.n_subjs,1);
for s=1:logg.n_subjs
    ws.d= subjdata{s,request.choicesession+1};
    ws.d = ws.d(ws.d(:, col.TrialOK)==1,:);
    ws.d(:, col.Subject)=s; 
    ws.d(:, col.PLdraw)= ws.d(:, col.PL(1)) ==ws.d(:, col.PL(2)); 
    ws.d(:, col.PPdraw)= ws.d(:, col.PP(1)) ==ws.d(:, col.PP(2)); 
   
% 	disp('Artificially use cf as stated - PILOT DATA ONLY');  ws.d(:, col.Trialcf)=0;  ws.d((ws.d(:, col.PL1) >ws.d(:, col.PL2)) & (ws.d(:, col.PP1) <ws.d(:, col.PP2)), col.Trialcf)=1; ws.d((ws.d(:, col.PL1) <ws.d(:, col.PL2)) & (ws.d(:, col.PP1) >ws.d(:, col.PP2)), col.Trialcf)=1;
    
%     
% 
%     % PP = residuals 
%     if s==1,  input('PP are residuals of PL. Not yet checked which variables this has been applied to.  Continue? '); end; 
%     [ws.b ws.bi ws.r]=  regress(ws.d(:, col.PP1),ws.d(:, col.PL1));
%     ws.d(:, col.PP1)= round( f_rescale(1, 10, ws.r));
%     [ws.b ws.bi ws.r]=  regress(ws.d(:, col.PP2),ws.d(:, col.PL2));
%     ws.d(:, col.PP2)=  round(f_rescale(1, 10, ws.r));


    % RT processing
    ws.rtout_high = (ws.d(:, col.ChoiceRT) > mean(ws.d(:, col.ChoiceRT)) + 2*std(ws.d(:, col.ChoiceRT))   )   ;
    ws.rtout_low = (ws.d(:, col.ChoiceRT) < mean(ws.d(:, col.ChoiceRT)) - 2*std(ws.d(:, col.ChoiceRT)) )   ;
    ws.d( (ws.rtout_high+ws.rtout_low)>0 , col.ChoiceRT)  =nan;
    ws.d(find(1- (ws.rtout_high+ws.rtout_low)), col.ChoiceRT) = frt(ws.d(find(1- (ws.rtout_high+ws.rtout_low)), col.ChoiceRT) );   % Preproc RTs 
    d_rt{s}=ws.d(:, col.ChoiceRT); 
    
    % Conflict - options are flagged     
    ws.d(:, col.PLcf) = fcf(  ws.d(:, col.PL(1)), ws.d(:, col.PL(2)) );
    ws.d(:, col.PPcf) = fcf(  ws.d(:, col.PP(1)), ws.d(:, col.PP(2)) );
    ws.d(:, col.PLPPcf) = fcf(  ws.d(:, col.PLPP(1)), ws.d(:, col.PLPP(2)) );
    
    subjdata{logg.n_subjs+1,request.choicesession+1}=[subjdata{logg.n_subjs+1,request.choicesession+1}; ws.d];    
    subjdata{s,request.choicesession+1}=ws.d;
    
    
    
    % Special samples
    d_ncf{s,request.choicesession+1} = ws.d(ws.d(:, col.Trialcf)==0,:);
    d_ncf{logg.n_subjs+1,request.choicesession+1}=[d_ncf{logg.n_subjs+1,request.choicesession+1};  d_ncf{s,request.choicesession+1}  ];    
    %
    ws.dncf=ws.d(ws.d(:, col.Trialcf)==0,:);
    d_ncf_nd{s,request.choicesession+1}= ws.dncf(ws.dncf(:, col.PLdraw)==0 & ws.dncf(:, col.PPdraw)==0, :);
    %     d_ncf_nd{s,request.choicesession+1} = ws.d(ws.d(:, col.Trialcf)==0 & ws.d(:, col.PLbest)~=0 & ws.d(:, col.PPbest)~=0 ,:);
    d_ncf_nd{logg.n_subjs+1,request.choicesession+1}=[d_ncf_nd{logg.n_subjs+1,request.choicesession+1};  d_ncf_nd{s,request.choicesession+1}  ];    
    %
    d_cf{s,request.choicesession+1} = ws.d(ws.d(:, col.Trialcf)==1,:);
    d_cf{logg.n_subjs+1,request.choicesession+1}=[d_cf{logg.n_subjs+1,request.choicesession+1};  d_cf{s,request.choicesession+1}  ];    
    %
    d_ncf_drawPL{s,request.choicesession+1}  =ws.d(ws.d(:, col.PLbest)==0, :); 
    d_ncf_drawPL{logg.n_subjs+1,request.choicesession+1}=[d_ncf_drawPL{logg.n_subjs+1,request.choicesession+1};  d_ncf_drawPL{s,request.choicesession+1}  ];    
    %
    d_ncf_drawPP{s,request.choicesession+1}  =ws.d(ws.d(:, col.PPbest)==0, :); 
    d_ncf_drawPP{logg.n_subjs+1,request.choicesession+1}=[d_ncf_drawPP{logg.n_subjs+1,request.choicesession+1};  d_ncf_drawPP{s,request.choicesession+1}  ];    
    %
    d_nd{s,request.choicesession+1} = ws.d(ws.d(:, col.PLbest)~=0 & ws.d(:, col.PPbest)~=0 ,:);
    d_nd{logg.n_subjs+1,request.choicesession+1}=[d_nd{logg.n_subjs+1,request.choicesession+1};  d_nd{s,request.choicesession+1}  ];    
    %
    d_dr{s,request.choicesession+1}= ws.d(ws.d(:, col.PLdraw)==1 | ws.d(:, col.PPdraw)==1 ,:);
    d_dr{logg.n_subjs+1,request.choicesession+1}=  [d_dr{logg.n_subjs+1,request.choicesession+1};  d_dr{s,request.choicesession+1}];
     
    
%     ws.d = sortrows(ws.d, [col.Trialcf col.PLbest col.PPbest]);
%     ws.d(:, [col.Trialcf col.PLbest col.PPbest])


    
    
    %
    ws=[]; 
end

% % [Request: Special samples ONLY?] ##########
% subjdata= d_ncf;   input('Requested NO conflict trials (draws not deleted) only!!! Continue?  ');
% subjdata= d_ncf_nd;   input('Requested NO conflict trials (no draws) only!!! Continue?  ');
% subjdata= d_ncf_drawPL;   input('Requested draw PL (no conflic) trials only!!! Continue?  ');
% subjdata= d_ncf_drawPP;   input('Requested draw PP (no conflic) trials only!!! Continue?  ');
subjdata= d_cf;   input('Requested beh-conflict trials (no draws) only!!! Continue? ');
% subjdata= d_nd;   input('Draws deleted!!! Continue?  ');
for s=1:logg.n_subjs;  d_rt{s}= subjdata{s, request.choicesession+1}(:, col.ChoiceRT); end

 
% input('THIS SCRIPT NEEEDS TO BE BUG CHECKED!'); 
sd = subjdata{1,3}



sd(:, [col.d1m2.PL col.d1m2.PP col.PLcf col.PPcf]) 
col

%% GLM settings 

request.ivs={  % Must match col contents  
    
    % [ Within-attribute comparison ] --------
%     'PL1'; 'PL2';'PP1'; 'PP2';     % Havent thorught through whether one needs to flip this for the RT analysis yet. Currently these are flipped
        'd1m2.PL'; 'd1m2.PP';     % choice 1>2, difference in each attribute
        'PLcf';  'PPcf';  
        
    % [ Linearly combining PL & PP scores ] --------
%     'PLPP1';  'PLPP2'; 
%     'd1m2.PLPP';     % choice 1>2 
%         'PLPPcf'; 

    % [ Event ratings ] --------
%     'er.d1m2.social'; 
%     'er.d1m2.happiness'; 
%     'er.d1m2.arousal';
%     'er.d1m2.experience';
%     'er.d1m2.vivid';
%     'er.d1m2.valence';
%     'er.ad1m2.social';  % absdiff
%     'er.ad1m2.happiness';
%     'er.ad1m2.vivid';
%     'er.ad1m2.valence';
%     'er.ad1m2.experience';
%     'er.ad1m2.arousal';
};



% Assemble IV columns + subject regressors (categorical)
col.ivs= cell2mat(cellfun(@(x,col)eval(['col.' x]), request.ivs, repmat({col}, length(request.ivs),1), 'UniformOutput',0)  );   % :D :D :D :D :D :D :D 
col.SubRegs_start = size(subjdata{logg.n_subjs+1,request.choicesession+1},2)+1;    
request.iv_names=cellfun(@(x)request.alliv_names(strcmp(request.alliv_names(:,1), x),2), request.ivs);   % if it crashes here, go updated the ist of iv names. DO NOT change this to a cell output 
subjdata{logg.n_subjs+1,request.choicesession+1}  = [subjdata{logg.n_subjs+1,request.choicesession+1} zeros(size(subjdata{logg.n_subjs+1,request.choicesession+1},1), logg.n_subjs)];
for s=1:logg.n_subjs ;  subjdata{logg.n_subjs+1,request.choicesession+1}(   subjdata{logg.n_subjs+1,request.choicesession+1}(:, col.Subject)==s,  col.SubRegs_start +s-1)=1;  end 
if request.IVs_rescale   % Re-scale IVs if requested
    error('i dont think you should do this. pl and pp should be relatable to each other')
    for s=1:logg.n_subjs+1
    ws.d =  subjdata{s, request.choicesession+1};  
    ws.d(:, col.ivs )=  (  ws.d(:, col.ivs ) - repmat(min(ws.d(:, col.ivs)), size(ws.d,1),1)  ) ./ repmat(max(  ws.d(:, col.ivs ) - repmat(min(ws.d(:, col.ivs)), size(ws.d,1),1)  ) , size(ws.d,1),1) ;  % Adjust floor, scale by ceiling 
    subjdata{s, request.choicesession+1}= ws.d ;
    end
end 
disp(' ------------------------------------------------------------------------------------------------------------------------')


%% Fix

do_fixed= 1; 
request.psig=0.08; 

if do_fixed
    % IVs 
    d_fixed=subjdata{logg.n_subjs+1,request.choicesession+1};  d_ivsfx = d_fixed(:, col.ivs); 
    ws.fx_subregs = d_fixed(:, col.SubRegs_start:col.SubRegs_start+logg.n_subjs-1) ;  ws.n_subregs = logg.n_subjs;   % Dummy subject regressors 
%         ws.fx_subregs = d_fixed(:, col.Subject);    ws.n_subregs = 1;   % Single subject regressor 
    disp('### Fixed effects results ##########')
    
    % Choice
    [b d st] = glmfit([ws.fx_subregs d_ivsfx], d_fixed(:, col.Choice)==1, 'binomial', 'link', 'logit');
    %     r_choicefx= [[{' ' 'Cst' } request.ivs']; [{'b'}  num2cell(b')]; [{'p'}  num2cell(st.p')]]; % sub regs will throw an error here 
    b(st.p> request.psig)=nan;  st.p(st.p> request.psig)=nan;
    disp('[Choice] IV         beta         pval -------- '); disp([request.ivs num2cell(b(ws.n_subregs +2:end)) num2cell(st.p(ws.n_subregs +2:end))])
    
    % RT       -  Non-conflict IVs must re-compiled (because difference scores are all Choice 1 >2)
    d_ivsfx_rt = d_ivsfx;  d_ivsfx_rt(d_fixed(:, col.Choice)~=1, find(cellfun(@(x)isempty(x), strfind(request.ivs, 'cf')) )) = -1 * d_ivsfx(d_fixed(:, col.Choice)~=1,find(cellfun(@(x)isempty(x), strfind(request.ivs, 'cf')) ));
    [b d st] = glmfit( [ws.fx_subregs  d_ivsfx_rt], d_fixed(:, col.ChoiceRT)    ); 
    %     r_rtfx= [[{' ' 'Cst'} request.ivs']; [{'b'}  num2cell(b')]; [{'p'}  num2cell(st.p')]]; % sub regs will throw an error here 
    b(st.p> request.psig)=nan;  st.p(st.p> request.psig)=nan;
    disp('[RT]  IV              beta         pval -------- '); disp([request.ivs num2cell(b(ws.n_subregs +2:end)) num2cell(st.p(ws.n_subregs +2:end))])
end


% error('done w fixedfx'); 

%% Mixed effects (1st level, 2nd level)

request.psig=0.08; 

b_Choicemx=nan(logg.n_subjs, length(col.ivs));  b_RTmx=b_Choicemx;
for s=1:logg.n_subjs  % First level
    ws.ivs = subjdata{s,request.choicesession+1}(:, col.ivs); ws.d= subjdata{s,request.choicesession+1};  
    
    % Choice 
    [ws.b ws.dd ws.st] = glmfit( ws.ivs , subjdata{s,request.choicesession+1}(:, col.Choice)==1, 'binomial', 'link', 'logit');
    b_Choicemx(s,:)= ws.b(2:end)';
    
    % RT 
    ws.ivs_rt=ws.ivs;  ws.ivs_rt(ws.d(:, col.Choice)~=1,  find(cellfun(@(x)isempty(x), strfind(request.ivs, 'cf')) ) ) = -1 *ws.ivs_rt(ws.d(:, col.Choice)~=1,  find(cellfun(@(x)isempty(x), strfind(request.ivs, 'cf')) ) );
    [ws.b ws.dd ws.st] = glmfit( ws.ivs_rt, ws.d(:, col.ChoiceRT) );
    b_RTmx(s,:)= ws.b(2:end)';
    ws=[]; 
end


% a=0
% openvar a
% 
% openvar b_Choicemx
% mean(b_Choicemx)
[h p ci st]=ttest(b_Choicemx); p

for i=1:size(b_Choicemx,2) 
    [p h st]=signrank(b_Choicemx(:,i));
    disp([request.iv_names{i} '  : p=' num2str(p)])
end
 

% 
% b_Choicemx-a
% 
% openvar ans

% Second levels + report 
request.dvs={'Choice';'RT'};
r_mx =  [[{' '} request.ivs']; [request.dvs cell(length(request.dvs), length(request.ivs))]]; 
for d=1:length(request.dvs)
    eval(['wd.betas = b_' request.dvs{d} 'mx;'])
    [wd.h wd.p wd.ci wd.st]=ttest(wd.betas  );
    for i=1:length(request.ivs)
        wi.st= num2str(wd.st.tstat(i),3);
        
        if wd.p(i)<0.001, wi.st=['***  p=' num2str(wd.p(i),3) ];
        elseif wd.p(i)<0.01, wi.st=['**  p=' num2str(wd.p(i),3) ];
        elseif wd.p(i)<0.05, wi.st=['*   p=' num2str(wd.p(i),3) ];
        elseif wd.p(i)<0.1, wi.st=['t   p=' num2str(wd.p(i),3) ];
        else wi.st=' ';
%             wi.st= ['    p=' num2str(wd.p(i),3) ];
        end
        r_mx{d+1, i+1}=wi.st;
    end
end  
disp('### Mixed effects results ##########'), r=r_mx'; disp(r_mx)
openvar r_mx

% Plot
do_fig=1;   do_plotindiv=1;
if do_fig, close all hidden
    f.plotcols=2;  f.maxplot=  2;
    f.figwidth= 800; f.figheight=400; f.fontsize=20; f.fontsize_title=30;f.fontname='PT Sans Caption';
    f.subplot_VerHorz=[0.15 0.1]; f.fig_BotTop=[0.15 0.15]; f.fig_LeftRight=[0.05 0.05];
    figure('Name', 'GLM betas', 'Position', [100 550 f.figwidth f.figheight], 'Color', 'w');  k=1;
    
    for d=1:length(request.dvs)
        eval(['wd.d= b_' request.dvs{d} 'mx;'])
        
        % Choice
        subtightplot(ceil(f.maxplot/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
        barwitherr(std(wd.d)./sqrt(logg.n_subjs),mean(wd.d),'y' );
        if do_plotindiv; hold on;  scatter(sortrows(repmat(1:length(col.ivs),1, logg.n_subjs)'),wd.d(:)); end 
        set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1:length(col.ivs) , 'XTickLabel', request.iv_names)
        title(request.dvs{d} ,  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
%         xticklabel_rotate, set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname)
        ylim([min(wd.d(:))  - std(wd.d(:)) max(wd.d(:))  + std(wd.d(:))])
        
    xlim([0.5 4.5])
        

%         % IV1>IV2
%         subtightplot(ceil( f.maxplot/ f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
%         barwitherr(std( wd.d(:,1) - wd.d(:,2) )./sqrt(logg.n_subjs),mean(wd.d(:,1) - wd.d(:,2)),'y' );
%         if do_plotindiv; hold on;  scatter(ones(logg.n_subjs,1),  wd.d(:,1) - wd.d(:,2)); end 
%         set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'TickDir','out',  'XTick', 1 , 'XTickLabel', [request.ivs{1} ' > ' request.ivs{2}] ,  'FontSize',f.fontsize, 'FontName', f.fontname);
%         title(request.dvs{d} ,  'FontWeight','bold', 'FontSize',f.fontsize_title, 'FontName', f.fontname,'Interpreter', 'none')
        
    end
end





% Plot rts (check transform distribution 
do_plotrt=0;
if do_plotrt
    f.plotcols=5;    f.figwidth= 500; f.figheight=1800; f.fontsize=10; f.fontsize_title=20;f.fontname='PT Sans Caption';
    f.subplot_VerHorz=[0.1 0.03]; f.fig_BotTop=[0.05 0.05]; f.fig_LeftRight=[0.05 0.05];
    figure('Name', 'Subject RT DVs', 'Position', [100 100 f.figwidth f.figheight], 'Color', 'w');  k=1;
    for s=1:logg.n_subjs
        subtightplot(ceil((logg.n_subjs+f.plotcols) / f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
        hist( d_rt{s}, 50 ); title(logg.subjects{s})
    end
    subtightplot(ceil((logg.n_subjs+f.plotcols) / f.plotcols),  f.plotcols, k:k+f.plotcols-1 ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);
    hist( d_rt{s}, 50 );  title(['All subs: ' logg.define{strcmp(logg.define(:,1), 'RT preproc') & cellfun(@(x)~isempty(x), logg.define(:,3)),2} ], 'FontSize',f.fontsize_title, 'FontName', f.fontname)
    if (strcmp(logg.define(:,1), 'RT preproc') & strcmp(logg.define(:,2), 'Log-transformed within subject'));  set(gca, 'xscale', 'log'); xlabel('RT (Log scale)');  end
    set(gca, 'FontSize',f.fontsize, 'FontName', f.fontname, 'box','on','xcolor',[0 0 0],'ycolor',[0 0 0])
end



%% End


disp('########################'), disp('Preproc decisions:');  disp(' '), disp([char(logg.define(:,1)) char(repmat({'  '}, size(logg.define,1), 1)) char(logg.define(:,2))] )

% max(d_ivsfx)
% min(d_ivsfx)

 
