% Plot model predictions (for comparison with data): Generate & plot choice contingencies PREDICTED by (specified) models
% clear all; close all hidden; clc; request.plotbins=[];  request.subsample=[];
clear all;  clc; request.plotbins=[];  request.subsample=[];

% Request 
request.plotbins={[1 2 3 4]; [5  6  7];  [8 9 10]};  % artificially bin simulated scores 
details.resfile='res_fit_s2_CC (17-Nov-2015)';   request.subsample='cc';
% details.resfile='res_fit_s2_IC (17-Nov-2015)';    request.subsample='ic';
% details.resfile='res_fit_s2_CC (19-Nov-2015) cf_decorraw'; 
% details.resfile='res_fit_s2_b3_CC (17-Nov-2015)';   request.plotbins=[]; 
request.n_simreps=500;

% Which models  
 for o=1:1 
    details.modfams.b=  {'b01' ; 'b02_l' };
    details.modfams.bc=  {'bc01' ; 'bc02_l' };
    details.modfams.heuristic={'h01_b_choPL'}; 
    
    % Append lapse families 
    details.modfams.bp= cellfun(@(x)[x(1) 'p' x(2:end)], details.modfams.b, 'UniformOutput',0); 
%     details.modfams.bpk= cellfun(@(x)[x(1) 'pk' x(2:end)], details.modfams.b, 'UniformOutput',0); 
 end
 
whichmodels=[1 ];  
details.whichmodels=[
%     details.modfams.b(whichmodels)
%     details.modfams.bp(whichmodels)
    details.modfams.bc(whichmodels)
%     details.modfams.heuristic(1);
    ];
 

for o1=1:1 % General settings and specifications
   
    % Paths 
    request.specificsubjects=[]; 
    where.where='/Users/EleanorL/Dropbox/SCRIPPS/7 Pleasure purpose';   
    where.mod =[where.where filesep '4b Modelling']; where.beh=[where.where filesep '3 Behaviour']; 
    path(pathdef), addpath(where.where), addpath([where.where fs '4a Beh analysis basic'])
    addpath([where.mod fs '1 Value functions' ])
    addpath([where.mod fs '1 Value functions' fs 'Heuristic'])
    
    % Model settings + fetch requested model details
%     request.results_folder 
    d_fits=load([where.mod filesep '2 Inputs' filesep  details.resfile ]);  details.fixedpar=d_fits.details.fixedpar;  details.n_iter= d_fits.details.n_iter; details.subj_ntrialsok= d_fits.details.subj_ntrialsok; 
    [details.model_defaults  details.par_transformations details.models] = f_modelsettings(details.whichmodels,details.n_iter);
    details.n_models=length(details.whichmodels);  details.fixedpar=[];  errorlog={}; e=1;
    w.options=optimset('Display', 'off', 'LargeScale','off'); rand('state',sum(100*clock));  % For fminunc
    details.dataset=d_fits.details.dataset;
    
    
    % Identify fit details + load data w subject selection 
    if isempty(strfind(details.resfile, 'pilot'))==0; details.session=1; else   details.session=2; end
    if isempty(strfind(details.resfile, '_b'))==0;   details.data_bin=num2str(details.resfile(strfind(details.resfile, '_b')+2));  else  details.data_bin=0; end
    if isempty(strfind(details.resfile, 'IC'))==0;    details.data_subset='IC'; elseif isempty(strfind(details.resfile, 'CC'))==0;  details.data_subset='CC';  else details.data_subset=[]; end
    w= load(['2 Inputs' fs '2 Data' fs details.dataset]);  % Manually assign dataset here if you want 
    col=w.logg.col;  details.define =w.logg.define;  details.dataset=w.logg.dataset; details.fxns.fcf= details.define{strcmp(details.define(:,1), 'Conflict'),3};
    eval(['col=col.s' num2str(details.session) ';']);  % Correct session columns
    %
    subjdata =w.subjdata;    details.subjects=subjdata(:,1); details.n_subjs=length(details.subjects);
    if isempty(request.specificsubjects)==0; error('not set up subject selection yet!'); end 
    details.subjects=subjdata(:,1); details.n_subjs=size(subjdata,1); 
    
    
    % Model settings + fetch requested model details (in requested order)
    d_fits.details4fit=d_fits.details; d_fits=rmfield(d_fits, 'details');
    orig.r_res=d_fits.r_res; 
    orig.r_iterations=[]; % disp('r_iterations turned off! all other lines w r_iterations also turned off. search r_iterations')
    orig.r_iterations=d_fits.r_iterations; % change to grid if needed
    if isempty(details.whichmodels); details.whichmodels=orig.r_res(1:5,1); end
    details.n_models=length(details.whichmodels); 
    d_fits.r_res=cell(length(details.whichmodels),  size(d_fits.r_res,2));
    d_fits.r_iterations=cell(length(details.whichmodels),  size(d_fits.r_iterations,2));
    for m=1:  length(details.whichmodels)
        if sum(strcmp(orig.r_res(:,1), details.whichmodels{m}))==1
            wm.rownum=find(strcmp(orig.r_res(:,1), details.whichmodels{m}));
            d_fits.r_res(m, 1:length(orig.r_res(wm.rownum,:)))= orig.r_res(wm.rownum,:);
            %
            wm.rownum=find(strcmp(orig.r_iterations(:,1), details.whichmodels{m}));
            d_fits.r_iterations(m, 1:length(orig.r_iterations(wm.rownum,:)))= orig.r_iterations(wm.rownum,:);
            %
            wm=[];
        else error(['Cannot find requested model: '         details.whichmodels{m}])
        end
    end
    
    % Interface
    [details.model_defaults  details.par_transformations details.models] = f_modelsettings(details.whichmodels);
    disp('====================================')
    disp(['Requested: ' num2str(size(subjdata,1)) ' subjects']); disp(' ')
    disp(['Fit from session : ' num2str(details.session)]);
    disp(['Fit res file: ' details.resfile]);  disp(' ')
    disp([ num2str(length(details.whichmodels)) ' models requested:']); disp(details.whichmodels(:)); disp(' ')
    % input('Hit enter to start                                   ');
    disp('====================================')
end

d_fits.r_res=sortrows(d_fits.r_res,3);
m=1;      L=  -1*sum(d_fits.r_res{m,2}(:,2) );   R =  sum(details.subj_ntrialsok).*log(1/2); disp(['Pseudo r2=  ' num2str( 1-(L/R) ) ' ,   p(Choice)=  '  num2str(  exp(L/sum(details.subj_ntrialsok)) )])
% Pseudo R2= 1-L/R, R=log L under chance. exp(L/ntrials) gives probability (relatable to chance)

%% (1) Generated choices as predicted by models 
%   d_fits: results from model fit

% Data variables 
disp(['Winning model:    '  d_fits.r_res{1,1} ]);
d_fits.original_res=d_fits.r_res;

% Artificially alter details
% disp('Parameters altered!!')
% d_fits.r_res{1,2}=d_fits. r_res{1,2}(1,:); details.subjects=details.subjects(1); details.n_subjs=1; disp('One subject only!');
% d_fits.r_res{1,2}(1,3+1:end)=[1 0 1 -12 0 0];   % cF 
% d_fits.r_res{1,2}(1,3+1:end)=[1 0 1 0 0 0];     % ct
% p=6;  d_fits.r_res{1,2}(:,3+p)=-0.1;
% d_fits.r_res{1,2}(:,4:end)=repmat([1 0 1   -12 0 0], 20,1);    % cF 
% d_fits.r_res{1,2}(:,4:end)=repmat([10 0 1 0  0 0], 20,1);
% d_fits.r_res{1,2}(:,4:end)=repmat([10 0 1 10], 1,1);

%% (1) Generate value & choices predicted by models 
frt=@(x) log(x);  
% frt=@(x) zscore(x);    


disp('Outstanding issue: intermittent discrepancies in L that I havent yet accounted for');

for o1=1:1 % Archive code 
% request.plotpar_steps=5;
% request.plotparrange={  % par num, par range
%     3   logspace(-1, 1 , request.plotpar_steps);    % bpma:  m param
% %     4   linspace(-5,15 , request.plotpar_steps);        % bpma: a param
%     4   linspace(-3,10 , request.plotpar_steps);        % bpma: a param
%     };
% 
% % compars(param, parval, model, choice, e, n)
% compars=nan(size(request.plotparrange,1),  request.plotpar_st eps, details.n_models, 3, 6,6);  % Variable for saving
% 
% 
% for iPar=1:size(request.plotparrange,1)
%     parnum=request.plotparrange{iPar,1};
%     parvals=request.plotparrange{iPar,2};
%     
%     for iVal=1:length(parvals)
%         d_fits.r_res=d_fits.original_res;
%         d_fits.r_res{1,2}(:,3+parnum)=parvals(iVal);  disp('Artificially altering parameters!!');
end


% Columns (for data to be added)
for o2=1:1
    
    % Existing (from plot_explorebeh)
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
%     col.cfweight=col.colmax+26;      % Ratio of PL:PP cf
%     col.cfweightbin=col.colmax+27;
%     col.cfweightlog=col.colmax+28;
%     col.cfweightbinrecip=col.colmax+29;
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
    col.PLcho=col.choPL;  col.PPcho=col.choPP;  col.PLuncho=col.unchoPL;  col.PPuncho=col.unchoPP;  col.PLPPcho=col.choPLPP;    col.PLPPuncho=col.unchoPLPP;
    col.PLchobin=col.choPLbin;  col.PPchobin=col.choPPbin;  col.PLunchobin=col.unchoPLbin;  col.PPunchobin=col.unchoPPbin;  col.PLPPchobin=col.choPLPPbin;    col.PLPPunchobin=col.unchoPLPPbin;
    
    % --------------------------------------------------------
    
    % Quantities that are specific to simulations 
    col.maxcol  = structmax(col);
    col.vOpt1 = col.maxcol  +1;
    col.vOpt2 = col.maxcol  +2;
    col.pOpt1 = col.maxcol  +3;
    col.pOpt2 = col.maxcol  +4;
    col.ObservedChoice= col.maxcol  +5;
end

% Generate simulated subjdata: d_pred{model, subject}
for s=1:details.n_subjs  % Alternations to subjdata itself
    ws.d=subjdata{s,details.session+1};
    ws.d=ws.d(ws.d(:, col.TrialOK)==1,:);  
    ws.d(:, col.ObservedChoice)=ws.d(:, col.Choice);
    
%     if s==1,  input('Regress raw scores out of cf?  ');  end; if isempty(strfind(details.resfile, ' cf_decorraw'))==1; input('ATTENTION! resfit file is not marked as cf-decorr-raw. Continue? '); end
%     [ws.b ws.bi ws.r]=  regress(ws.d(:, col.PLcf),ws.d(:, [col.PL(1) col.PL(2)])); ws.d(:, col.PLcf)=ws.r;    
%     [ws.b ws.bi ws.r]=  regress(ws.d(:, col.PPcf),ws.d(:, [col.PP(1) col.PP(2)])); ws.d(:, col.PPcf)=ws.r;
%     
    % RT pre-processing
    ws.rtout_high = (ws.d(:, col.ChoiceRT) > mean(ws.d(:, col.ChoiceRT)) + 2*std(ws.d(:, col.ChoiceRT))   )   ;
    ws.rtout_low = (ws.d(:, col.ChoiceRT) < mean(ws.d(:, col.ChoiceRT)) - 2*std(ws.d(:, col.ChoiceRT)) )   ;
    ws.d( (ws.rtout_high+ws.rtout_low)>0 , col.ChoiceRT)  =nan;
    ws.d(find(1- (ws.rtout_high+ws.rtout_low)), col.ChoiceRT) = frt( ws.d(find(1- (ws.rtout_high+ws.rtout_low)), col.ChoiceRT) );   % Preproc RTs
    
    % Random calcs  
    
    
    subjdata{s,details.session+1}=ws.d;
end


for m= 1:details.n_models  % Implement simulations
    disp(['Model ' num2str(m) ':  '   details.models{m} ' ---------------------------'])
    d_pred{details.n_subjs+1,m}=[];
    wm.model=details.whichmodels{m};
    wm.modrow_res=find(strcmp(d_fits.r_res(:,1), details.whichmodels{m}));
    %             disp(['Sum of nLLs: ' num2str(sum(d_fits.r_res{wm.modrow_res,2}(:,2)))])
    
    % Subject-specific predictions followed by mean predictions
    for s=1:details.n_subjs
        ws.modpar=d_fits.r_res{wm.modrow_res,2}(s,4:end);  ws.origpar=ws.modpar; % Fetch subject's model parameters
        ws.d=subjdata{s,details.session+1}; 
        
        % Apply inverse transforms
        for p=1:details.models{m,2}
            x=ws.modpar(p); eval(['ws.modpar(p)= ' details.models{m, 5}{p}  ';'])
        end
        if isempty(strfind(details.models{m,1}, 'p'))==1;   [nll pch]=f_nllsoftmax(ws.modpar, {wm.model ws.d details.fixedpar col}); 
        else [nll pch]=f_nllsoftmax_lapse(ws.modpar, {details.models{m,1} ws.d details.fixedpar col});
            % [nll pch]=f_nllsoftmax_lapse(ws.modpar, {details.models{m,1} ws.d details.fixedpar col});
        end;  ws.d(:, col.Choice)=nan; 
        if  d_fits.r_res{m,2}(s,2)-nll~=0;   disp(['nLL calc discrepancy!    '  details.models{m,1} '  (' subjdata{s,1} ')    '  num2str(d_fits.r_res{m,2}(s,2),4) ' vs ' num2str(nll, 4)]); end % Paranoid check
        
        % V(Choices) according to this model
        eval(['ws.v=' wm.model '(ws.modpar,  {[] ws.d details.fixedpar col});'])
        ws.d(:, [col.vOpt1  col.vOpt2])=squeeze(ws.v) ;
        
        % Calculate p(Choice) using softmax rule
        ws.beta=ws.origpar(1);  ws.softmaxbase=( exp(ws.beta.*ws.d(:, col.vOpt1))+exp(ws.beta.*ws.d(:, col.vOpt2) ));
        if isempty(strfind(details.models{m,1}, 'p'))~=1; disp('havent worked out how to implement softmax choice for p!')
            ws.epsilon=ws.origpar(2); ws.softmaxbase= ws.softmaxbase .*(1-2.*ws.epsilon + ws.epsilon);
        end
        ws.d(:,col.pOpt1)=exp(ws.beta.*ws.d(:, col.vOpt1)) ./ ws.softmaxbase;
        ws.d(:,col.pOpt2)=exp(ws.beta.*ws.d(:, col.vOpt2)) ./ ws.softmaxbase;
        
        % Implement predicted choie + calculate subsequent variables
        if s==details.n_subjs; disp('  (Choice generated acc to models)'); end 
        ws.d=repmat(ws.d, request.n_simreps, 1);  ws.randcho= rand(size(ws.d,1),1);  ws.d(ws.d(:,col.pOpt1)<= ws.randcho , col.Choice)=1;  ws.d(ws.d(:,col.pOpt1)>ws.randcho , col.Choice)=2;  % Stochastic prediction         
        ws.d= f_calcfromraw(ws.d, col);
        ws.d=f_calcnew(ws.d, col, request.plotbins, details.fxns);
        
        for o2=1:1 % Hard-coded simulations (checks)
%         % % [ Hard-coded simulations ] ------
%         disp('Hard code acc to conflict'); ws.dd=ws.d; ws.dd(:, col.Choice)=nan; 
%         if unique(ws.d(:, col.PLcf))>3;  ws.dd(:,  [col.PLcf col.PPcf] )= ws.dd(:,  [col.PLbincf col.PPbincf]);  end % For raw/unbinned scores. Assume cf values of 3-5
%         ws.wh=ws.dd(:, col.PLcf)==5 & ws.dd(:, col.PPcf)==5;  ws.dd(ws.wh, col.Choice)=randi(3, sum(ws.wh),1);
%         ws.wh = (ws.dd(:, col.PLcf)==4 & ws.dd(:, col.PPcf)==5)+ (ws.dd(:, col.PLcf)==5 & ws.dd(:, col.PPcf)==4);  ws.wh =find(ws.wh);
%             ws.dd(ws.wh(1:round(length(ws.wh)*0.6)), col.Choice)=  ws.dd(ws.wh(1:round(length(ws.wh)*0.6)), col.PLPPbest);
%             ws.dd(ws.wh(round(length(ws.wh)*0.6)+1:end), col.Choice)=  abs(ws.dd(ws.wh(round(length(ws.wh)*0.6)+1:end), col.PLPPbest)-3);
%         ws.wh= (ws.dd(:, col.PLcf)==3 & ws.dd(:, col.PPcf)==5)+ (ws.dd(:, col.PLcf)==5 & ws.dd(:, col.PPcf)==3) + (ws.dd(:, col.PLcf)==4 & ws.dd(:, col.PPcf)==4);  ws.wh =find(ws.wh);
%             ws.dd(ws.wh(1:round(length(ws.wh)*0.7)), col.Choice)=  ws.dd(ws.wh(1:round(length(ws.wh)*0.7)), col.PLPPbest);
%             ws.dd(ws.wh(round(length(ws.wh)*0.7)+1:end), col.Choice)= abs(ws.dd(ws.wh(round(length(ws.wh)*0.7)+1:end), col.PLPPbest)-3);
%         ws.wh= (ws.dd(:, col.PLcf)==3 & ws.dd(:, col.PPcf)==4)+ (ws.dd(:, col.PLcf)==4 & ws.dd(:, col.PPcf)==3);  ws.wh =find(ws.wh);
%             ws.dd(ws.wh(1:round(length(ws.wh)*0.8)), col.Choice)=  ws.dd(ws.wh(1:round(length(ws.wh)*0.8)), col.PLPPbest);
%             ws.dd(ws.wh(round(length(ws.wh)*0.8)+1:end), col.Choice)= abs(ws.dd(ws.wh(round(length(ws.wh)*0.8)+1:end), col.PLPPbest)-3);
%         ws.wh=ws.dd(:, col.PLcf)==3 & ws.dd(:, col.PPcf)==3;  ws.dd(ws.wh, col.Choice)= ws.dd(ws.wh, col.PLPPbest);
%         ws.d(:, col.Choice)=ws.dd(:, col.Choice);
%         ws.d= f_calcfromraw(ws.d, col);
%         ws.d=f_calcnew(ws.d, col, request.plotbins, details.fxns);
        end
        
        
        % Carry through choice
        ws.d(:, [col.RTchoPL col.RTchoPP col.RTchoPLPP])=nan;
        ws.d(ws.d(:, col.Trialcf)==1 & ws.d(:, col.conchoPL)==1,  col.RTchoPL)= ws.d(ws.d(:, col.Trialcf)==1 & ws.d(:, col.conchoPL)==1,  col.ChoiceRT);
        ws.d(ws.d(:, col.Trialcf)==1 & ws.d(:, col.conchoPP)==1,  col.RTchoPP)= ws.d(ws.d(:, col.Trialcf)==1 & ws.d(:, col.conchoPP)==1,  col.ChoiceRT);
        ws.d(ws.d(:, col.Trialcf)==1 & ws.d(:, col.conchoPLPP)==1,  col.RTchoPLPP)= ws.d(ws.d(:, col.Trialcf)==1 & ws.d(:, col.conchoPLPP)==1,  col.ChoiceRT);

        
        %  HUGE redundancy here in terms of empty matrices. Kept none the less so that different models can be plotted at the same time.
        %       Data that you want to plot needs to be allocated here
        
        
        
        
        % ######################################################################
        % INSERT HERE To use subject's model parameters to calculate other variables (theoretical best? Rank order choices?)
        
        
        
        % ######################################################################
        d_pred{s,m}=ws.d;
        d_pred{details.n_subjs+1,m}=[d_pred{details.n_subjs+1,m}; d_pred{s,m}];
        ws=[];
    end
    
    %
    wm=[];
end


for o1=1:1 % Archive code
    
    %         % Record values for params : compars(param, parval, model, choice, e, n)
    %         for m=1:details.n_models
    %             for c=1:3
    %                 compars(iPar, iVal, m, c, 1:6, 1:6)= squeeze( v_choicematrix{details.n_subjs+1, m, c});
    %             end
    %         end
    %     end
    % end
    
    
    % Plot par ranges of compars(param, parval, model, choice, e, n)
    % for m=1:details.n_models
    %     for iPar=1:size(request.plotparrange,1);
    %         figure('Name',[ '[' details.models{m, 1}      '] Choice/vals over range param ' details.models{m,3}{request.plotparrange{iPar,1}}], 'NumberTitle', 'off', 'Position', [130 85 900 600]);         set(gcf,'Color',[1 1 1])
    %         parnum=request.plotparrange{iPar,1};
    %         parvals=request.plotparrange{iPar,2};
    %
    %         for iVal=1:length(parvals);
    %             for c=1:3
    %                 subplot(request.plotpar_steps,3, (iVal-1)*3+c);
    %                 imagesc(squeeze(compars(iPar, iVal, m, c, :, :))); axis square; colorbar;
    %
    %
    % %                 imagesc(squeeze(compars(1, iVal, m, c, :, :))+squeeze(compars(2, iVal, m, c, :, :)) ); axis square; colorbar;
    %
    %
    %                 caxis([ min(min(min(squeeze(compars(iPar, :, m, c, :, :)))))            max(max(max(squeeze(compars(iPar, :, m, c, :, :))))) ]);
    % %                 caxis([ min(min(min(min(squeeze(compars(:, :, m, c, :, :))))))            max(max(max(max(squeeze(compars(:, :, m, c, :, :)))))) ]);
    %             end
    %         end
    %     end
    % end
    
    % % Difference?
    % figure;
    %
    % steps=[1 0 ]; k=1;
    % steps=[3 1 ]; k=1;
    % for ss=1:2
    %
    % subplot(4,2, k); imagesc(squeeze(compars(1, steps(1)+ss, 1, 3, :, :)   -  compars(1, steps(2)+ss, 1, 3, :, :))); axis square; axis off; k=k+1;
    %
    % title([ details.models{m,3}{request.plotparrange{1,1}} ' param: step '  num2str(steps(1)+ss) '  - step ' num2str(steps(2)+ss)])
    % subplot(4,2,k); imagesc(squeeze(compars(2, steps(1)+ss, 1, 3, :, :)   - compars(2, steps(2)+ss, 1, 3, :, :))); axis square; axis off;  k=k+1;
    % title([ details.models{m,3}{request.plotparrange{2,1}} ' param: step '  num2str(steps(1)+ss) ' - step ' num2str(steps(2)+ss)])
    % end
end

% Which subsample? ------
% request.subsample='all';
% request.subsample='cc';
% request.subsample='ic';
% request.subsample='cc_nd';
% request.subsample='nd';
for o=1:1
    input(['[Subsample] Requested: ' request.subsample '.  Continue?  ']); 
    switch request.subsample
        case 'cc'
            d_pred=cellfun(@(x, col)x(x(:,col)==0,:),  d_pred,repmat({col.Trialcf}, size(d_pred)),'UniformOutput',0);
            if strcmp(details.data_subset,'CC')==0; disp('Requesting CC sub-sample on fit NOT run on CC'); end 
        case 'cc_nd'
            d_pred=cellfun(@(x, col)x(  x(:,col(1))==0 & x(:,col(2))==0 & x(:,col(3))==0,:),  d_pred,repmat({[col.Trialcf col.PLdraw col.PPdraw]}, size(d_pred)),'UniformOutput',0);
            if strcmp(details.data_subset,'CC')==0; disp('Requesting CC sub-sample on fit NOT run on CC'); end 
        case 'ic'
            d_pred=cellfun(@(x, col)x(x(:,col)==1,:),  d_pred,repmat({col.Trialcf}, size(d_pred)),'UniformOutput',0);
            if strcmp(details.data_subset,'IC')==0; disp('Requesting IC sub-sample on fit NOT run on IC'); end 
        case 'nd'
            d_pred=cellfun(@(x, col)x(x(:,col(1))==0 & x(:,col(2))==0,:),  d_pred,repmat({[ col.PLdraw col.PPdraw]}, size(d_pred)),'UniformOutput',0);
    end
end

% Insert here: sub-sample trials



%% [2] Plot requested x for PLdimension vs PP dimension
% Always PL first, PP second 

% Request
request.ivs={
    'bincf';
    %         'marbin';
    'chobin';
    };
request.ivs={'bincf'; 'chobin';  'marbin'; };  % Binned 
% request.ivs={'cf'; 'cho';  'mar'; };          

request.dvs =  {  % DVs: If any of these need to be calculated on the fly (i.e. within specified IV cells), make sure they are computed below!
    'pChoPLPP'; 
%     'pChoPL'; 'pChoPP';
    %         'pChoPLPPbin'; 'pChoPLbin'; 'pChoPPbin';
    %         'ChoiceRT'; 'RTchoPL'; 'RTchoPP'
    } ;

% Calculate scores for requested binning 
d_subd=repmat({cell(2,2)}, [details.n_models,details.n_subjs, size(request.ivs,1)]);   % {m,subject, iv}{PLbin, PPbin}: all data for that cell (for subject)
d_subntrials=repmat({nan(2,2)}, details.n_subjs+1, size(request.ivs,1));   % {subject, iv}{PLbin, PPbin}: n trials for that cell (for subject)
d_dv=repmat({cell(2,2, details.n_subjs )}, [details.n_models, size(request.ivs,1),   size(request.dvs,1)]);   % {m,iv, dv}{PLbin, PPbin, subject}: all data for that cell (for subject)
for m=1:details.n_models  % Calculate scores 
 
    for o=1:1  % Setup
        % [Columns] IVs are assumed to have PL & PP dimensions; DVs are assumed to not.
        request.iv_cols{m} =  cell2mat( cellfun(@(x, col)[eval(['col.PL' x]) eval(['col.PP' x])], request.ivs(:,1) ,  repmat({col},  size(request.ivs,1),1),  'UniformOutput',0) );
        request.dv_cols{m} =   cell2mat( cellfun(@(x, col)eval(['col.' x]), request.dvs,  repmat({col},size(request.dvs,1),1), 'UniformOutput',0));
        
        % IV bins
        for iv=1:size(request.ivs,1)
            request.ivs{iv,3}= {unique( d_pred{details.n_subjs+1,m}(:, request.iv_cols{m}(iv,1)) ) unique( d_pred{details.n_subjs+1,m}(:, request.iv_cols{m}(iv,2)) )};
            request.ivs{iv,3}{1} = request.ivs{iv,3}{1}(find(1-isnan(request.ivs{iv,3}{1})));
            request.ivs{iv,3}{2}=request.ivs{iv,3}{2}(find(1-isnan(request.ivs{iv,3}{2})));
            if m==1, d_subntrials{details.n_subjs+1, iv}= zeros(length(request.ivs{iv,3}{1}), length(request.ivs{iv,3}{2})); end
        end
        
    end
    
    for s=1:details.n_subjs % Extract data
        ws.d= d_pred{s,m};
        %     ws.d=sortrows(ws.d, [col.PPcf col.PLcf]);
        
        
        % Bin raw data into IV cells + calculate new quantities if you wan
        for iv=1:size(request.ivs,1)
            for ipl=1:length(request.ivs{iv,3}{1})
                for ipp=1:length(request.ivs{iv,3}{2})
                    wc.d= ws.d(   ws.d(:,  request.iv_cols{m}(iv, 1)) == request.ivs{iv,3}{1}(ipl) & ws.d(:,  request.iv_cols{m}(iv, 2)) == request.ivs{iv,3}{2}(ipp) , :);
                    if m==1, d_subntrials{s,iv}(ipl, ipp)=size(wc.d,1); d_subntrials{details.n_subjs+1,iv}(ipl, ipp)=d_subntrials{details.n_subjs+1,iv}(ipl, ipp)+  d_subntrials{s,iv}(ipl, ipp); end
                    %                     wc.d(:, col.nTrials)=size(wc.d,1);
                    
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
                        
                        
                        
                        d_subd{m,s,iv}{ipl, ipp} = wc.d;
                    else d_subd{m,s,iv}{ipl, ipp} = nan;  % If you change this, it will add zeros to all unfilled columns
                    end
                    
                end
            end
        end
        
        % Auto-read DVs
        for iv=1:size(request.ivs,1)
            for ipl=1:length(request.ivs{iv,3}{1})
                for ipp=1:length(request.ivs{iv,3}{2})
                    for dv=1:size(request.dvs,1)  % Auto-read DVs for non-empty cells (no on-the-fly calculation
                        if isempty(d_subd{m,s,iv}{ipl, ipp} )==1 | (isnan(d_subd{m,s,iv}{ipl, ipp}) & size(d_subd{m,s,iv}{ipl, ipp},2)==1)   % empty d_subd cells had a nan inserted
                            d_dv{m,iv,dv}{ipl, ipp, s} =nan;
                        else d_dv{m,iv,dv}{ipl, ipp, s} = nanmean(d_subd{m,s,iv}{ipl, ipp}(:, request.dv_cols{m}(dv)));
                        end
                    end
                end
            end
        end
        
        ws=[];
    end
    
end 

% Plot
%     f.axis=[1.5 1.9];   % Comment out to omit
% close all , 
f.plotcols=  size(request.dvs,1); % f.plotcols= 4;
f.figwidth= 800; f.figheight=500; f.fontsize=20; f.fontsize_title=20;f.fontname='PT Sans Caption';
    f.subplot_VerHorz=[0.15 0.05]; f.fig_BotTop=[0.15 0.1]; f.fig_LeftRight=[0.2 0.2];
% f.subplot_VerHorz=[0.15 0.015]; f.fig_BotTop=[0.05 0.05]; f.fig_LeftRight=[0.25 0.25];
for m=1:details.n_models
    figure('Name',details.models{m}, 'Position', [100+m*30 250-m*30 f.figwidth f.figheight], 'Color', 'w'); k=1;
    for iv=1:size(request.ivs,1)
        for dv=1:size(request.dvs,1)
            wc.d = d_dv{m,iv,dv}; wc.plot = nanmean(cell2mat(wc.d) ,3); % This is what you want 
            
            if isempty(strfind(request.dvs{dv}, 'RT'))==0; wc.nancol=[0 0.4 0]; else wc.nancol=[0 0.0 0]; end
            subtightplot(size(request.ivs,1), f.plotcols , k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
            imagescnan(wc.plot  , 'nancolor', wc.nancol), colorbar, axis square
%             title(['[' request.ivs{iv} ']   ' request.dvs{dv}],'FontSize', f.fontsize)
            title( request.dvs{dv},'FontSize', f.fontsize), ylabel(request.ivs{iv,1},'FontSize', f.fontsize_title)
            set(gca,'FontSize', f.fontsize)
%                     ylabel(['PL ' request.ivs{iv}],'FontSize', f.fontsize), xlabel(['PP ' request.ivs{iv}],'FontSize', f.fontsize)
            set(gca, 'ytick',  1:length(request.ivs{iv,3}{1}), 'yticklabel', num2cell(request.ivs{iv,3}{1}),  'xtick',  1:length(request.ivs{iv,3}{2}),  'xticklabel', num2cell(request.ivs{iv,3}{2}))
            
            if isempty(strfind(request.dvs{dv}, 'RT'))==0; colormap 'pink'; end
            if isfield(f,'axis'),  caxis(f.axis),    disp('Artificial cxis!'),end
            
             % Mark which model
            if iv==1 && dv==1; text(-1.6,4.5, ['Model: ' details.models{m}], 'Fontsize',f.fontsize_title,'rotation',90); end
           
        end
    end
end


% Plot individual
do_plotindiv= 0;
if do_plotindiv
    iv=1; dv=1;       f.plotcols= 5; f.figwidth= 800; f.figheight=250; f.fontsize=13; f.fontsize_title=13;f.fontname='PT Sans Caption';
    f.subplot_VerHorz=[0.05 0.01]; f.fig_BotTop=[0.05 0.05]; f.fig_LeftRight=[0.15 0.25];   
    for m=1:details.n_models
         figure('Name', ['Individual sims:  ' details.models{m}], 'Position', [800 0 f.figwidth f.figheight], 'Color', 'w'); k=1;
        for s=1:details.n_subjs
           
            subtightplot(ceil(details.n_subjs/f.plotcols) , f.plotcols , k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1;
            imagescnan(cell2mat(d_dv{m,iv,dv}(:,:,s)) , 'nancolor', [0 0 0]), colorbar, axis square
%                 title(details.subjects{s},'FontSize', f.fontsize)
            axis off, set(gca,'FontSize', f.fontsize)
            %     ylabel(['PL ' request.ivs{iv}],'FontSize', f.fontsize), xlabel(['PP ' request.ivs{iv}],'FontSize', f.fontsize)
            
            
            % Mark which model
            if s==1;   text(-0.2,15.4, ['[ ' details.models{m} ' ]  IV: ' request.ivs{iv,1} '   DV: ' request.dvs{dv}], 'Fontsize',f.fontsize_title,'rotation',90);   end
        end
    end
end
    
    
%% End


