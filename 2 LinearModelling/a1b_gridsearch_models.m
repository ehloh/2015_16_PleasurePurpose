% Fit models to the choice data, calculate BICs & parameter values
clear all;close all hidden; clc

% Request specific
details.n_iter=10; 
details.session=2;    % 2= fMRI, 3= in lab 
%
details.specificsubjects={'p01_YH'}; % BLANK to process all subjects
% details.specificsubjects='AllDataOK';  % Only complete data ok
details.dataset='All data (11-Nov-2015)';

% Which models  
 for o=1:1 
    details.modfams.b=  {'b01' ; 'b02_l' };
    details.modfams.bc=  {'bc01' ; 'bc02_l' };
    details.modfams.heuristic={'h01_b_choPL'}; 
    
    % Append lapse families 
    details.modfams.bp= cellfun(@(x)[x(1) 'p' x(2:end)], details.modfams.b, 'UniformOutput',0); 
%     details.modfams.bpk= cellfun(@(x)[x(1) 'pk' x(2:end)], details.modfams.b, 'UniformOutput',0); 
    
    
end


details.whichmodels=details.modfams.bp(1);
details.whichmodels=details.modfams.heuristic(1); 
details.whichmodels=[
%     details.modfams.b; 
    details.modfams.bc; 
%     details.modfams.bp
    ];

for o1=1:1 % General settings and specifications
   
    % Paths 
    where.where='/Users/EleanorL/Dropbox/SCRIPPS/7 Pleasure purpose';   
    where.mod =[where.where filesep '4b Modelling']; where.beh=[where.where filesep '3 Behaviour']; 
    path(pathdef), addpath(where.where)
    addpath([where.mod fs '1 Value functions' ])
    addpath([where.mod fs '1 Value functions' fs 'Heuristic'])
    
    % Load subjects
    [n t r]=xlsread([where.beh filesep 'datalog_plpr.xlsx']); 
    w= load(['2 Inputs' fs '2 Data' fs details.dataset]);  subjdata =w.subjdata;  col=w.logg.col; details.define =w.logg.define;  details.dataset=w.logg.dataset; 
    details.subjects=subjdata(:,1); details.n_subjs=length(details.subjects);
    if isempty(strfind(details.dataset, 'pilot'))==1 & isempty(strfind(details.dataset, 'Pilot'))==1
        [n t r]=xlsread([where.beh filesep 'datalog_plpr.xlsx']);
        [logg.subjects logg.n_subjs logg.datalog] = f_selectsubjects(r, details.specificsubjects,r, 'All');
    else
        logg.subjects=[cellfun(@(x)['p0' num2str(x)], num2cell(1:9),'UniformOutput',0)'; cellfun(@(x)['p' num2str(x)], num2cell(10:20),'UniformOutput',0)';];
        logg.n_subjs = length(logg.subjects); 
        [logg.subjects logg.n_subjs logg.datalog] = f_selectsubjects([[{'Subject'}; logg.subjects]  [{'All'}; num2cell(ones(20,1))] ], details.specificsubjects,[[{'Subject'}; logg.subjects]  [{'All'}; num2cell(ones(20,1))] ], 'All');
    end
     
    % Model settings + fetch requested model details
    [details.model_defaults  details.par_transformations details.models] = f_modelsettings(details.whichmodels,details.n_iter);
    details.n_models=length(details.whichmodels);  details.fixedpar=[];  errorlog={}; e=1;
    w.options=optimset('Display', 'off', 'LargeScale','off'); rand('state',sum(100*clock));  % For fminunc
    diary([where.mod fs '2 Inputs' filesep '1 Fit logs' filesep 'diary_fit_' num2str(details.session) ' (' date ')'])
    eval(['col=col.s' num2str(details.session) ';']);  % Correct session columns
    details.subj_ntrials = cellfun(@(x)size(x,1), subjdata(:, details.session+1 ) );
    details.subj_ntrialsok = cellfun(@(x)size(x,1),   cellfun(@(x)x(x(:, col.TrialOK)==1, :),  subjdata(:, details.session+1 ),'UniformOutput',0 ));
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
    disp(['No. of subjects: ' num2str(details.n_subjs)])
%     if isempty(details.specificsubjects)==0;  if details.n_subjs == 19 & sum(strcmp(details.subjects, {'p01_YH';'p02_MI';'p03_AY';'p04_AW';'p05_CA';'p06_BT';'p07_HC';'p08_KC';'p09_KJ';'p10_YC';'p11_BI';'p13_MS';'p14_SK';'p15_MM';'p16_SH';'p17_BB';'p18_WB';'p19_HB';'p20_LZ';}) ) ==19; disp('   Subjects w complete data only (n=19)'); else disp('   Subset of subjects only:'); disp(details.specificsubjects); end;   end; disp(' ')
    disp(['Requested: ' num2str(details.n_subjs) ' subjects, ' num2str(details.n_iter) ' steps per param'])
    disp(['Choice session:  ' details.session]); disp(' ');
    disp([num2str(length(details.whichmodels))  ' Models:']); disp(details.whichmodels); disp(' ');
    disp('Paramater ranges:'); for p=1:size(details.par_transformations,1); disp(['      ' details.par_transformations{p,1} ':      ' num2str(details.par_transformations{p,4}(1)) '     to    ' num2str(details.par_transformations{p,4}(end))]); end; disp(' ');
    disp(' '); input('Hit Enter to start      ')
    disp('=======================================================')
    
end 

%% Grid search all parameter combinations



startm=1;

% Fit grid: 'r_grid' - model grids for each subject
%       Col 1=Model name
%       Col 2=Grid search results (nLL, parameters)
%             Col 1: Subject
%             Col 2: nLL
%             Col 3: Full grid of nLL (all param-points, p-dimensions)
%             Col 4: Indices of best fit (calculated later)
%             Col 6 onwards: best fit model parameters (1st parameter is beta/inverse temperature)
if exist('r_grid', 'var')==0; r_grid=cell(details.n_models,1); end
for m=startm:details.n_models
    disp(['Model ' num2str(m) ' - ' details.models{m,1} '  ##############' ])
    r_grid{m,1}=details.models{m,1};
        
    % Apply inverse constraints to parameter ranges (correct for constraints
    %       being applied within the softmax/value functions). Col 7=Parameters
    %       with inverse constraints
    wi.startpar=cell(1,details.models{m,2});
    for p=1:details.models{m,2}
        x=details.models{m,6}{p};
        eval(['wi.startpar{p}=' details.models{m,5}{p} ';']);
    end
    details.models{m,7}=wi.startpar;
    
    % Fit all subjects
    for s=1: details.n_subjs
        disp(['Subject ' num2str(s) '  (' details.subjects{s} ')'])
        r_grid{m,2}{s,1}=details.subjects{s};
        
        if strcmp(details.models{m,1}(2:3), 'pk')==1; 
            eval(['[ ws.nLL, ws.np, ws.ns ] = f_optimgrid_vect(@f_nllsoftmax_sklapse, {'''  details.models{m,1} ''' subjdata{s, details.session+1}( subjdata{s, details.session+1}(:,col.TrialOK)==1, :)  details.fixedpar col},  details.models{m,7});'])
        elseif strcmp(details.models{m,1}(2), 'p')==1; eval(['[ ws.nLL, ws.np, ws.ns ] = f_optimgrid_vect(@f_nllsoftmax_lapse, {'''  details.models{m,1} ''' subjdata{s, details.session+1}( subjdata{s, details.session+1}(:,col.TrialOK)==1, :)  details.fixedpar col},  details.models{m,7});'])
        else eval(['[ ws.nLL, ws.np, ws.ns ] = f_optimgrid_vect(@f_nllsoftmax, {'''  details.models{m,1} ''' subjdata{s, details.session+1}( subjdata{s, details.session+1}(:,col.TrialOK)==1, :)  details.fixedpar col},  details.models{m,7});'])            
%              [ws.nLL, ws.np, ws.ns ] = f_optimgrid_vect(@f_nllsoftmax, {'b01' subjdata{s, details.session+1}( subjdata{s, details.session+1}(:,col.TrialOK)==1, :)  details.fixedpar col},  details.models{m,7});
        end
        r_grid{m,2}{s,3}=ws.nLL;
    end
    
end


%% For each subject, find best fit (+ calculate model vals)

% 'r_res'
%       Col 1: Model name
%       Col 2: Subject fit parameters
%             Col i: BIC
%             Col ii: nLL
%             Col ii: Pseudo r2
%             Col iv onwards: parameters (beta first)
%       Col 3: Model BIC (summed across subjects)
%       Col 4:
%       Col 5:
%       Col 6 onwards: mean parameter values
r_res=cell(details.n_models,  5+max(cell2mat(details.models(1,2))));
for o1=1:1 % Columns 
    rc.modname=1;
    rc.subpars=2;
    rc.sp.bic=1;
    rc.sp.nll=2;
    rc.sp.p1=4;
    rc.modelbic=3;
    rc.mean_p1=6;
    %
    r_res=cell(details.n_models,  5+max(cell2mat(details.models(1,2)))); r_iterations=cell(details.n_models,1); 
end
for m= 1:details.n_models
    disp(['Model ' num2str(m) ' - ' details.models{m,1} '  ##############' ])
    r_res{m,1}=details.models{m,1};
    r_res{m,2}=nan(details.n_subjs, 3+details.models{m,2});
    
    wm.indices=[]; % Indices for best fit
    for i=1: details.models{m,2}; wm.indices= [wm.indices ' ii' num2str(i)]; end
    
    for s=1:details.n_subjs
        disp(['Subject ' num2str(s) '  (' details.subjects{s} ')'])
        r_grid{m,2}{s,1}=details.subjects{s};
        ws.gridnLL=r_grid{m,2}{s,3};
        
        % Locate best match
        ws.nmatch=length(find(ws.gridnLL==min(ws.gridnLL(:))));
        eval(['[' wm.indices ']=ind2sub(size(ws.gridnLL), find(ws.gridnLL==min(ws.gridnLL(:))));'])
        if ws.nmatch~=1; disp('Multiple matches! Assume 1st (arbitrary)'); for i=1: details.models{m,2}; eval(['ii' num2str(i) '=ii' num2str(i) '(1);']); end; end
        eval(['ws.index=num2cell([' wm.indices ']);']);
        ws.ii=  [cellfun(@(x)[num2str(x) ','], ws.index(1:length(ws.index)-1), 'UniformOutput',0) num2str(ws.index{length(ws.index)}) ];
        r_grid{m,2}{s,4}=ws.index;
        
        % Record details of best fit (r_grid)
        eval(['r_grid{m,2}{s,2}=ws.gridnLL(' strcat(ws.ii{:}) ');'])
        for p=1:details.models{m,2}
            r_grid{m,2}{s,5+p}=details.models{m,6}{p}(ws.index{p});
        end
        
        % Record details of best fit (r_res{m,2})
        r_res{m,2}(s,2)=r_grid{m,2}{s,2};
        r_res{m,2}(s,3)= 1 -  (-1*r_res{m,2}(s,2) ./ details.subj_ntrialsok(s).*log(0.5));  %Pseudo R2= 1-L/R, R=log L under chance
        
        r_res{m,2}(s,4:3+details.models{m,2})=cell2mat(r_grid{m,2}(s,6:end));
        r_res{m,2}(s,1)= 2*r_res{m,2}(s,2) +  details.models{m,2} * log(size(subjdata{s,details.session+1},1)); % Calculate BIC: 2*nll+ K*ln(n_trials)
    end
    
    % Calculate overall BIC & mean parameters for all subjects
    r_res{m,3}=sum(r_res{m,2}(:,1));
    r_res(m, 6:5+details.models{m,2})=num2cell(mean(r_res{m,2}(:,4:end),1));
    
end
r_res=sortrows(r_res,3);

% Save
resfilename=['res_gridmodels_s' num2str(details.session) ' (' date ')'];
resfilewhere=[where.mod filesep '2 Inputs' filesep ];
resfilesthere=cellstr(spm_select('List', resfilewhere, ['res_gridmodels_s' num2str(details.session) '*.*' date]));
if isempty(resfilesthere)~=1 && strcmp(resfilesthere{end}(strfind(resfilesthere{end}, ')')+1), '.')==0 % number re existing fit from this date
    k=1; kk=0;
    while kk==0 
        if sum(strcmp(resfilesthere,[resfilename num2str(k) '.mat']));  k=k+1;
        else resfilename=[resfilename num2str(k)]; kk=1;
        end
    end
end
save([resfilewhere resfilename],  'r_grid','r_res', 'details')

%% END 

disp('===================================='); w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
disp('INSTRUCTION: Results variables (prefix r):r_grid, r_res, r_models')
disp('                     [r_res: Col 1= Model, Col 2=Subject fits, 3=Model BIC, Col 4=N valid subjects, Col 6 onwards=Mean parameter values]'); disp(' ')
disp('Check command-window log for poor fits; See ''further instructions'' at end of script for Bayes factor & model weights'); disp(' ')
disp('====================================')


%% Further instructions (Copy to command window to calculate

% Calculate Bayes factor ----------------------
m1=1;  m2=2;
B=(r_res{m1,3}-r_res{m2,3})*-0.5;

disp(['B=' num2str(B)  '  (m1=' r_res{m1,1} ', m2='  r_res{m2,1} ')'])
openvar r_res;
% Interpreting B (conventions): 3-10=moderate evidence, >10=strong evidence (in favour of m1)

try f_sendemail('kurzlich', ['Modelling fitting (using grid) done [' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ']'], ' ',1); end
%
% % ----------------------------------------------


% Fits pars?
m=1; disp(['Model:   '   r_grid{m,1} '     ' details.models{m,1}     '       '  r_res{m,1} ' (if mismatch-DELETE!)'])
parnum=2; pars=r_res{m,2}(:, 3+parnum);
hist(pars); details.models{m,6}{parnum}

%% Look directly at nLL surface
%       Chance nLL: nTrials x log(1/2)

m=1; s=15; 

if strcmp(r_grid{m,1}, details.models{m,1})+strcmp(r_grid{m,1}, r_res{m,1})~=2; error('Check order of models!'); end
disp(['Model name:   '  details.models{m,1}]);
gg=r_grid{m,2}; g=gg{s,3};
for p=4: size(r_res{m,2},2);  % Which grid point for best fit?
    disp([details.models{1,3}{p-3} ':  '  num2str(find(details.models{1,6}{p-3}==r_res{m,2}(s,p))) ])
end

% Plot
ggg =squeeze(g(10,:, :));
surf(ggg)
ylabel('Par 1'); xlabel('Par 2');  
title([r_grid{m,1} '  (' r_grid{m,2}{s,1} ', min='   num2str(r_res{m,2}(s,2), 4) ')'])


