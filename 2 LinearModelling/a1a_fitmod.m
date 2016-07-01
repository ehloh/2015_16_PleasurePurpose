% Fit models to choice data, calculate BIC and param values 
clear all;close all hidden; clc; details.res_suffix=[];

% Request specific
details.n_iter=20; 
details.session=2;    % 2= fMRI, 3= in lab 
%
details.specificsubjects={};  
 
% details.dataset='All data (11-Nov-2015)';
details.dataset='All data IC (17-Nov-2015)';
% details.dataset='All data CC (17-Nov-2015)';
% details.dataset='All pilot data (11-Nov-2015)'; details.session=1;    
% details.dataset='All pilot data IC (11-Nov-2015)'; details.session=1;    
% details.dataset='All pilot data CC (11-Nov-2015)'; details.session=1;    

% Which models  
 for o=1:1 
    details.modfams.b=  {'b01' ; 'b02_l' };
    details.modfams.bc=  {'bc01' ; 'bc02_l' };
    details.modfams.heuristic={'h01_b_choPL'}; 
    
    % Append lapse families 
    details.modfams.bp= cellfun(@(x)[x(1) 'p' x(2:end)], details.modfams.b, 'UniformOutput',0); 
    details.modfams.bpc= cellfun(@(x)[x(1) 'p' x(2:end)], details.modfams.bc, 'UniformOutput',0);  
 end

whichmodels=[1];   
details.whichmodels=details.modfams.bp(whichmodels);
details.whichmodels=details.modfams.heuristic(whichmodels); 
details.whichmodels=[
    details.modfams.b; 
    details.modfams.bc; 
%     details.modfams.bp;    
%     details.modfams.bpc;
    details.modfams.heuristic
    ];
 
for o1=1:1 % General settings and specifications
   
    % Paths 
    details.whichmodels= details.whichmodels(:);
    where.where='/Users/EleanorL/Dropbox/SCRIPPS/7 Pleasure purpose';   
    where.mod =[where.where filesep '4b Modelling']; where.beh=[where.where filesep '3 Behaviour']; 
    path(pathdef), addpath(where.where)
    addpath([where.mod fs '1 Value functions' ])
    addpath([where.mod fs '1 Value functions' fs  'Heuristic'])   
%     addpath([where.mod fs '1 Value functions' fs  'weight'])
%     
    % Load subjects
    [n t r]=xlsread([where.beh filesep 'datalog_plpr.xlsx']); 
    w= load(['2 Inputs' fs '2 Data' fs details.dataset '.mat']);   subjdata =w.subjdata;   if isfield(w, 'logg')==0; error('You forgot to run the prep modelling script!'), end;  
    col=w.logg.col;  details.define =w.logg.define;  details.subjects=subjdata(:,1); details.n_subjs=length(details.subjects);
    if isempty(strfind(details.dataset, 'pilot'))==1 & isempty(strfind(details.dataset, 'Pilot'))==1
        [n t r]=xlsread([where.beh filesep 'datalog_plpr.xlsx']);
        [logg.subjects logg.n_subjs logg.datalog] = f_selectsubjects(r, details.specificsubjects,r, 'All');
    else
        logg.subjects=[cellfun(@(x)['p0' num2str(x)], num2cell(1:9),'UniformOutput',0)'; cellfun(@(x)['p' num2str(x)], num2cell(10:20),'UniformOutput',0)';];
        logg.n_subjs = length(logg.subjects); 
        [logg.subjects logg.n_subjs logg.datalog] = f_selectsubjects([[{'Subject'}; logg.subjects]  [{'All'}; num2cell(ones(20,1))] ], details.specificsubjects,[[{'Subject'}; logg.subjects]  [{'All'}; num2cell(ones(20,1))] ], 'All');
    end
    if isempty(details.specificsubjects)==0; error('not set up subject selection yet!'); end 

    
    % Model settings + fetch requested model details
    [details.model_defaults  details.par_transformations details.models] = f_modelsettings(details.whichmodels,details.n_iter);
    details.n_models=length(details.whichmodels);  details.fixedpar=[];  errorlog={}; e=1;
    w.options=optimset('Display', 'off', 'LargeScale','off'); rand('state',sum(100*clock));  % For fminunc
    diary([where.mod fs '2 Inputs' filesep '1 Fit logs' filesep 'diary_fit_' num2str(details.session) ' (' date ')'])
    eval(['col=col.s' num2str(details.session) ';']);  % Correct session columns
    
    details.subj_ntrials = cellfun(@(x)size(x,1), subjdata(:, details.session+1 ) );
    details.subj_ntrialsok = cellfun(@(x)size(x,1),  cellfun(@(x)x(x(:, col.TrialOK)==1, :),  subjdata(:, details.session+1 ), 'UniformOutput',0 ) );
    startm=1;
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
    disp(['No. of subjects: ' num2str(details.n_subjs)])
    if isempty(details.specificsubjects)==0;  if details.n_subjs == 19 & sum(strcmp(details.subjects, {'p01_YH';'p02_MI';'p03_AY';'p04_AW';'p05_CA';'p06_BT';'p07_HC';'p08_KC';'p09_KJ';'p10_YC';'p11_BI';'p13_MS';'p14_SK';'p15_MM';'p16_SH';'p17_BB';'p18_WB';'p19_HB';'p20_LZ';}) ) ==19; disp('   Subjects w complete data only (n=19)'); else disp('   Subset of subjects only:'); disp(details.specificsubjects); end;   end; disp(' ')
    disp(['Requested: ' num2str(details.n_subjs) ' subjects, ' num2str(details.n_iter) ' iterations per model'])
    disp(['Choice session:  ' details.session]); disp(' ');
    disp([num2str(length(details.whichmodels))  ' Models:']); disp(details.whichmodels); disp(' ');
    disp('Paramater ranges:'); for p=1:size(details.par_transformations,1); disp(['      ' details.par_transformations{p,1} ':      ' num2str(details.par_transformations{p,4}(1)) '     to    ' num2str(details.par_transformations{p,4}(end))]); end; disp(' ');
    disp(' '); input('Hit Enter to start      ')
    disp('=======================================================')
    
end


%%


% Data alterations (marked in file suffix & details.define)
for s=1:details.n_subjs
    ws.d= subjdata{s,details.session+1};
        
%     if s==1,  input('Regress raw scores out of cf?  '); details.res_suffix=' cf_decorraw';  details.define(size(details.define,1)+1, 1:3)= {'cf scores have had opt1/2 scores regressed out'   'cf-decorr-raw' ' '}; end
%     [ws.b ws.bi ws.r]=  regress(ws.d(:, col.PLcf),ws.d(:, [col.PL(1) col.PL(2)])); ws.d(:, col.PLcf)=ws.r;    
%     [ws.b ws.bi ws.r]=  regress(ws.d(:, col.PPcf),ws.d(:, [col.PP(1) col.PP(2)])); ws.d(:, col.PPcf)=ws.r;
%     
    subjdata{s,details.session}=ws.d;
end


% input('Subjects all, turn on trycatch, turn on save!!!');
% startm=1;

% 'r_res'
%       Col 1: Model name
%       Col 2: Subject fit parameters
%             Col i: BIC
%             Col ii: nLL
%             Col iii: fminunc exceeded default iterations?
%             Col iv onwards: parameters (beta first)
%       Col 3: Model BIC (summed across subjects)
%       Col 4: Hessians
%       Col 5: 
%       Col 6 onwards: mean parameter values
for o1=1:1 % Columns 
    rc.modname=1;
    rc.subpars=2;
    rc.sp.bic=1;
    rc.sp.nll=2;
    rc.sp.p1=4;
    rc.modelbic=3;
    rc.hessians=4;
    rc.mean_p1=6;
    %
    if exist('r_res', 'var')==0; if startm~=1; input('Requested start is NOT from model #1. Proceed?    '); end; r_res=cell(details.n_models,  5+max(cell2mat(details.models(1,2)))); r_iterations=cell(details.n_models,4);  end
end
for m=startm:details.n_models
    w.c=clock; disp(['Model ' num2str(m) ' - ' details.models{m,1} '       ['  num2str(w.c(4))  ':'  num2str(w.c(5))  '] ##############' ])
    r_res{m,1}=details.models{m,1};
    r_iterations{m,1}=details.models{m,1};
    r_iterations{m,2}=zeros(details.n_subjs*details.n_iter, details.models{m,2}+2);
    r_iterations{m,3}=cell(details.n_subjs*details.n_iter, 2);
    details.models{m, 6}= cellfun(@(x)x(randperm(details.n_iter) ), details.models{m, 6}, 'UniformOutput',0);   % Randomize order of par steps

    % 'r_iterations' - model fits for each subject x iteration
    %       Col 1=Model name
    %       Col 2=Iteration results (nLL, parameters)
    %             Col 1: Subject
    %             Col 2: nLL
	%             Col 3: Pseudor2 / Iterations exceeded fminunc?
    %             Col 4 onwards: model parameters (1st parameter is beta/inverse temperature)
    %       Col 3=Iteration hessians
    %             Col 2: Hessian for this iteration
    %       Col 4=nLL histogram data
    %       Col 5= Subject hessians for best fit (col 1), ok? (non-singular; col 2)
    for s= 1: size(subjdata,1)
        disp(['Subject ' num2str(s) '  (' details.subjects{s} ')'])
        ws.data= subjdata{s, details.session+1}(subjdata{s, details.session+1}(:, col.TrialOK)==1, :);
        
%         ff=inline('-log( (20./x)-1)')
%         ff(3)
         
        for i=1:details.n_iter
            try 
                %  Starting parameters (read seed + apply inverse constraint if applicable)
                wi.startpar=nan(1,details.models{m,2});
                for p=1:details.models{m,2}
                    x=details.models{m,6}{p}(i); eval(['wi.startpar(p)=' details.models{m,5}{p} ';']);
                    if strcmp(details.models{m,3}{p}, 'm')~=1;  wi.startpar(p)=round(wi.startpar(p)*10)/10; end
                end
                
                % Fit parameters ####################################
                if strcmp(details.models{m,1}(2), 'p')==1
                    [wi.par, wi.L, wi.exit, wi.output, wi.grad, wi.hessians]=fminunc(@(x)f_nllsoftmax_lapse(x, {details.models{m,1} ws.data details.fixedpar col},0,0), wi.startpar, w.options);
                else
                    [wi.par, wi.L, wi.exit, wi.output, wi.grad, wi.hessians]=fminunc(@(x)f_nllsoftmax(x, {details.models{m,1} ws.data details.fixedpar col},0,0), wi.startpar, w.options);
%                     [nll pch]=f_nllsoftmax(wi.startpar, {details.models{m,1} ws.data details.fixedpar col});
                end
                
                % Convert parameters to parameter space
                for p=1:details.models{m,2}
                    x=wi.par(p); eval(['wi.par(p)=' details.models{m,4}{p} ';']);
                end
                
                % Write to array
                wi.rownum=(s-1)*details.n_iter+i;
                r_iterations{m,2}(wi.rownum,1)=s;
                r_iterations{m,2}(wi.rownum,2)=wi.L;
                r_iterations{m,2}(wi.rownum,3)= wi.exit==0;
                r_iterations{m,2}(wi.rownum,4:3+length(wi.par))=wi.par;
                r_iterations{m,3}{wi.rownum,1}=s;
                r_iterations{m,3}{wi.rownum,2}=wi.hessians;
                wi=[];
                
                catch MExc
                    errorlog{e,1}=['Failed: ' details.subjects{s} '  -  ' details.models{m,1} '    iteration ' num2str(i)];
                    errorlog{e,2}=wi.startpar;
                    errorlog{e,3}=MExc;
                    disp(errorlog{e,1}); e=e+1;
                    
                    % Fake inputs
                    wi.rownum=(s-1)*details.n_iter+i;
                    r_iterations{m,2}(wi.rownum,1)=s;
                    r_iterations{m,2}(wi.rownum,2)=nan;
                    r_iterations{m,2}(wi.rownum,4:3+length(wi.startpar))=wi.startpar;
                    r_iterations{m,3}{wi.rownum,1}=s;
                    r_iterations{m,3}{wi.rownum,2}=[];
            end
        end
        
        % Choose best iteration for this subject + record details
        ws.iters=r_iterations{m,2}(  (s-1)*details.n_iter+1:s*details.n_iter, :);
        ws.iterhess=r_iterations{m,3}((s-1)*details.n_iter+1:s*details.n_iter, :);
        ws.bestiter= ws.iters(find(ws.iters(:,2)== min(ws.iters(:,2)), 1, 'first'),:);
        r_res{m,2}(s,2)=ws.bestiter(2);
        r_res{m,2}(s,4:3+length(ws.bestiter)-3)=ws.bestiter(4:end);
        if ws.bestiter(3)==1;  r_res{m,2}(s, 3)=1; end
        r_res{m,4}{s,1}=ws.iterhess{find(ws.iters(:,2)== min(ws.iters(:,2)), 1, 'first'),2};
        r_iterations{m,5}{s,1}=ws.iterhess{find(ws.iters(:,2)== min(ws.iters(:,2)), 1, 'first'),2};
        try r_iterations{m,5}{s,2}= sum(isnan(r_iterations{m,5}{s,1}(:)))>0 + rank(r_iterations{m,5}{s,1})<details.models{m,2}; catch; r_iterations{m,5}{s,2} =nan; end
            
        % Calculate BIC: 2*nll+ K*ln(n_trials)
        r_res{m,2}(s,1)= 2*r_res{m,2}(s,2) +  details.models{m,2} * log(size(ws.data,1));
        r_res{m,2}(s,3)= 1 -  (-1*r_res{m,2}(s,2) ./ (details.subj_ntrialsok(s).*log(0.5)));  %Pseudo R2= 1-L/R, R=log L under chance
        %
        ws=[];
    end
    
    % Calculate overall BIC & mean parameters for all subjects 
    r_res{m,3}=sum(r_res{m,2}(:,1));
    r_res(m, 6:5+details.models{m,2})=num2cell(mean(r_res{m,2}(:,4:end),1));
    r_iterations{m,4}= mean(reshape(r_iterations{m,2}(:,2), details.n_iter, details.n_subjs),2);  % nLL histograms
    
    % Save in partial fits in-between models
    wi.c=clock; save(['2 Inputs' filesep 'PartialFits' filesep 'Partial fit s' num2str(details.session) ' (' date ' - ' num2str(wi.c(4)) ' hrs)']);
end
r_res=sortrows(r_res,3);

% try  f_sendemail('kurzlich', ['[' ThisPC '] Runthru done [' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ']'], ' ',1); end; char(errorlog{:,1})
% error('Continue to save?')


%% End

disp('===================================='); w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
disp('INSTRUCTION: Results variables (prefix r):r_iterations, r_res, r_models')
disp('                     [r_res: Col 1= Model, Col 2=Subject fits, 3=Model BIC, Col 4=N valid subjects, Col 6 onwards=Mean parameter values]'); disp(' ')
disp('Check command-window log for poor fits; See ''further instructions'' at end of script for Bayes factor & model weights'); disp(' ')
disp('====================================')


% Save
if isempty(strfind(details.dataset, 'pilot'))==0; 
resfilename=['res_fit_pilot']; % if isempty(strfind(details.dataset, 'data b'))==0;   resfilename=[resfilename 'b'  details.dataset(strfind(details.dataset, 'data b')+6)];  end
else resfilename=['res_fit_s' num2str(details.session)]; 
end
if isempty(strfind(details.dataset, 'data b'))==0;   resfilename=[resfilename '_b'  details.dataset(strfind(details.dataset, 'data b')+6)];  end
if isempty(strfind(details.dataset, 'IC'))==0;  resfilename=[resfilename '_IC'];  elseif isempty(strfind(details.dataset, 'CC'))==0;  resfilename=[resfilename '_CC']; end
resfilename=[resfilename ' (' date ')' details.res_suffix];
resfilewhere=[where.mod filesep '2 Inputs' filesep ]; cd(resfilewhere), filelist=dir; filelist=cellstr(char(filelist(3:end).name));  % filelist=cellstr(ls);
if sum(strcmp(filelist, [resfilename '.mat']))>0;  k=2; kk=0;
    while kk==0, if sum(strcmp(filelist, [resfilename  num2str(k) '.mat']))==0, resfilename=[resfilename num2str(k)]; kk=1; else k=k+1; end; end
end
save([resfilewhere resfilename], 'details', 'r_iterations','r_res', 'rc','errorlog'); diary off
try % Transfer file to DeletedDaily 
    movetowhere='\\Asia\DeletedDaily\EL ModRuns';
    if isdir(movetowhere)==0; mkdir(movetowhere);end
    w.c=clock; cd(resfilewhere); 
    transferfilename=[ThisPC '_' num2str(w.c(3)) '_' num2str(w.c(2)) 'at' num2str(w.c(4)) num2str(w.c(5)) 'hrs ' resfilename];
    copyfile([resfilename '.mat'],  [movetowhere filesep transferfilename '.mat']);        
    disp(['Saved locally:  '  resfilename ]); disp(['Saved in DeletedDaily:  '  transferfilename])
catch
    disp(['Saved locally:  '  resfilename ' - NOT transferred to DeletedDaily']);
end
try  f_sendemail('kurzlich', ['[' ThisPC '] Modelling fitting (using fminunc) done [' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ']'], ' ',1); end; char(errorlog{:,1})



%% Further instructions (Copy to command window to calculate)

% nLL histograms: can we trust the minimizations? Iterations should reliably find the lowest nLL fit
f.subplot_rowcols=[4 4];
% f.subplot_rowcols=[4 2];
f.figwidth= 500; f.figheight=1000; ff=1; mm=1;
f.subplot_VerHorz=[0.1 0.03]; f.fig_BotTop=[0.05 0.03]; f.fig_LeftRight=[0.05 0.1];
figure('Name', ['nLL Histograms Plot ' num2str(ff)], 'NumberTitle', 'off', 'Position',[200,00,f.figheight,f.figwidth], 'Color',[1 1 1]);
for m=1:details.n_models
    subtightplot(f.subplot_rowcols(1), f.subplot_rowcols(2),   mm,  f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);
    hist(round(r_iterations{m,4})); title(r_iterations{m,1}, 'FontSize', 14)
%     axis off
    if m~=details.n_models && m/(f.subplot_rowcols(1)*f.subplot_rowcols(2)) == round(m/(f.subplot_rowcols(1)*f.subplot_rowcols(2)))
        ff=ff+1; mm=1; 
        figure('Name', ['nLL Histograms Plot ' num2str(ff)], 'NumberTitle', 'off', 'Position',[200,00,f.figheight,f.figheight], 'Color',[1 1 1]);
    else mm=mm+1; 
    end
end

% Instead of BICs, look at likelihoods 
Ls =  cellfun(@(x)sum(x(:,2)), r_res(:,2));
 

m1=1;
m2=2;

% Calculate Bayes factor ---------------------- 
%        Which models to compare? (Row number in 'r_res')
L=  -1*sum(r_res{m1,2}(:,2) );   R =  sum(details.subj_ntrialsok).*log(1/2);
disp(['Pseudo r2=  ' num2str( 1-(L/R) ) ' ,   p(Choice)=  '  num2str(  exp(L/sum(details.subj_ntrialsok)) )])
B=(r_res{m1,3}-r_res{m2,3})*-0.5;
disp(['B=' num2str(B)  '  (m1=' r_res{m1,1} ', m2='  r_res{m2,1} ')'])
% Interpreting B (conventions): 3-10=moderate evidence, >10=strong evidence (in favour of m1)


%% Results BICs

if sum(cellfun(@(x)sum(x(:,2)<0), r_res(:,2)))>0;  disp('Do models have problematic Ls? (0s=OK)'), disp([r_res(:,1) num2cell(cellfun(@(x)sum(x(:,2)<0), r_res(:,2)))]); end 

for o1=1:1  % Figure settings 
fontsize=15;
fontname='PT Sans Caption';  % pt serif (caption) ,san serif , pt sans,trebuchet ms
% fontname='Cambria';
% fontname='Arial';
end

% Plot BICs
figure('Position', [100 -250 1000 600], 'Color', 'w');  
bar(cell2mat(sortrows(r_res(:,3), -1)));
set(gca,'FontSize',fontsize, 'FontName', fontname, 'LineWidth', 0.8,'TickDir','out'); 
xlabel('Model', 'FontSize',fontsize, 'FontName', fontname); ylabel('BIC', 'FontSize',fontsize, 'FontName', fontname)
ylim([7000 17000]);  
title('BICs across model space', 'FontSize',fontsize,  'FontName', fontname); 
xlim([0 65])

% title('BICs across model space in Control task', 'FontSize',fontsize, 'FontName', fontname);  xlim([0 65])
if isempty(details.define)==0;  disp('########################'), disp('Preproc decisions:');  disp(' '), disp([char(details.define(:,1)) char(repmat({'  '}, size(details.define,1), 1)) char(details.define(:,2))] ); end

% Exclude certain subjects & re-calculate BICs
include=[1:15 17:20];   input('DO: Subset subjects?') ; 
for m=1:size(r_res,1)
    r_res{m,2}= r_res{m,2}(include,:);
    r_res{m,3}=sum(r_res{m,2}(:,1));
end
details.subj_ntrialsok= details.subj_ntrialsok(include); 
r_res=sortrows(r_res,3); 
