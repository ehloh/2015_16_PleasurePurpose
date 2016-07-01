% Hierarchical LCA  
clear all;close all hidden; clc;  
details.whichPC='mac';


% Request specific
details.n_parsteps = 8 ;   % For grid search 
details.session=2;    % 2= fMRI, 3= in lab 
details.normalize_withinattr=1;  % Normalize scores within attribute (on a trial-level basis)
%
details.specificsubjects={};  

details.specificsubjects= {
%     'p01_YH'; 
% 'p02_MI'; 'p03_AY';'p04_AW';'p05_CA';
%     'p06_BT';'p07_HC';'p08_KC';'p09_KJ';
'p10_YC';
%     'p11_BI';'p12_AL';'p13_MS';'p14_SK'; 'p15_MM';
%     'p16_SH';'p17_BB';'p18_WB';'p19_HB';'p20_LZ';
    };

 
% details.dataset='All data (17-Nov-2015)';
% details.dataset='All data IC (17-Nov-2015)';
details.dataset='All data CC (17-Nov-2015)';
% details.dataset='All pilot data (11-Nov-2015)'; details.session=1;    
% details.dataset='All pilot data IC (11-Nov-2015)'; details.session=1;    
% details.dataset='All pilot data CC (11-Nov-2015)'; details.session=1;    

for o1=1:1 % Request models 
    details.modelfams{1}={'dltan01'}; 
    details.modelfams{2}={'lta01'};  
end
details.whichmodels=[details.modelfams{1}; details.modelfams{2} ];
details.whichmodels ={'lta01'};  

for o1=1:1 % General settings and specifications
   
    % Paths 
    w=pwd; if strcmp(w(1),'/'), where.where='/Users/EleanorL/Dropbox/SCRIPPS/3 Pleasure purpose';    else  where.where='C:\Users\e.loh\Dropbox\SCRIPPS\3 Pleasure purpose'; end 
    where.mod =[where.where filesep '4b Modelling']; where.beh=[where.where filesep '3 Behaviour']; 
    path(pathdef), addpath(where.where), cd(where.mod)
    addpath([where.mod filesep '1 Value functions'])  
    addpath([where.mod filesep '1 Value functions' filesep 'hLCA'])  
    
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
    if isempty(details.specificsubjects)==0;  [details.subjects details.n_subjs subjdata] = f_selectsubjects([ cell(1, size(subjdata,2) );  subjdata], details.specificsubjects, [[{'Subject'}; subjdata(:,1)]  [{'All'}; num2cell(ones(length(subjdata(:,1)),1))] ], 'All');  subjdata=subjdata(2:end,:);  end 
    w.options=optimset('Display', 'off', 'LargeScale','off'); rand('state',sum(100*clock));  % For fminunc
   
    % Model settings + fetch requested model details
    [details.model_defaults  details.par_transformations details.models] = f_hLCAsettings(details.whichmodels,details.n_parsteps );
    details.n_models=length(details.whichmodels);  details.fixedpar=[];  errorlog={}; e=1;
    w.options=optimset('Display', 'off', 'LargeScale','off'); rand('state',sum(100*clock));  % For fminunc
    eval(['col=col.s' num2str(details.session) ';']);  % Correct session columns
    fmat= @(mat, index)mat(index);
    rng(0, 'twister'); % to control variance of random no. generation 
    if isempty(strfind(details.dataset,'IC')), details.subset_trials='IC';
    elseif isempty(strfind(details.dataset,'CC')), details.subset_trials='CC';
    else details.subset_trials=[];
    end

    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
    disp(['No. of subjects: ' num2str(details.n_subjs)])
    if isempty(details.specificsubjects)==0;  if details.n_subjs == 19 & sum(strcmp(details.subjects, {'p01_YH';'p02_MI';'p03_AY';'p04_AW';'p05_CA';'p06_BT';'p07_HC';'p08_KC';'p09_KJ';'p10_YC';'p11_BI';'p13_MS';'p14_SK';'p15_MM';'p16_SH';'p17_BB';'p18_WB';'p19_HB';'p20_LZ';}) ) ==19; disp('   Subjects w complete data only (n=19)'); else disp('   Subset of subjects only:'); disp(details.specificsubjects); end;   end; disp(' ')
    disp(['Requested: ' num2str(details.n_subjs) ' subjects  '])
    disp(['Choice session:  ' details.session]); disp(' ');
    disp(['Subset trials:  ' details.subset_trials]); disp(' ');
    disp([num2str(length(details.whichmodels))  ' Models:']); disp(details.whichmodels); disp(' ');
    disp('Paramater ranges:'); for p=1:size(details.par_transformations,1); disp(['      ' details.par_transformations{p,1} ':      ' num2str(details.par_transformations{p,4}(1)) '     to    ' num2str(details.par_transformations{p,4}(end))]); end; disp(' ');
    disp(' '); input('Hit Enter to start      ')
    disp('=======================================================')
    
end  

%% Calculate grids for each subject 

req.do_subgrids= 1;      startm=1; 
if req.do_subgrids
    for o=1:1  % Settings for fit 
        fp.maxtime= 7000;
        fp.timebin =10;  % Size of time bin (ms)
        fp.nsamples=1000;  % Disabled for now
        
        % Columns for fed-in data
        fp.col.pl1=1;
        fp.col.pp1=2;
        fp.col.pl2=3;
        fp.col.pp2=4;
        fp.col.ch=5;
        fp.col.rt=6;
        for p=1:length(details.model_defaults{strcmp(details.model_defaults(:,1), 'hLCA'),3})
            if ~isempty(details.par_transformations{strcmp(details.par_transformations(:,1), details.model_defaults{strcmp(details.model_defaults(:,1), 'hLCA'),3}{p}), 5}),  eval(['fp.'  details.model_defaults{strcmp(details.model_defaults(:,1), 'hLCA'),3}{p} ' =' num2str(details.par_transformations{strcmp(details.par_transformations(:,1), details.model_defaults{strcmp(details.model_defaults(:,1), 'hLCA'),3}{p}), 5}) ';']);  end
        end
        details.fixedpar=fp; 
        disp(['Requested grid: '  num2str(details.n_parsteps) ' param steps, ' num2str(fp.nsamples) ' samples per run'])
    end
    for s=1: details.n_subjs  % Set up data 
        ws.dd= subjdata{s, details.session+1};
        
        % Get rid of trials in which there are no ratings/choice 
        ws.dd = ws.dd(~isnan(sum(ws.dd(:, [ col.PL(1) col.PP(1) col.PL(2) col.PP(2)]),2)),:);  
        ws.dd= ws.dd(ws.dd(:, col.Choice)~=0, :);

        % Normalize within-attribute
        if details.normalize_withinattr 
            ws.dd(:, [ col.PL(1) col.PL(2)]) = ws.dd(:, [ col.PL(1) col.PL(2)])./repmat(sum(ws.dd(:, [ col.PL(1) col.PL(2)]),2),1,2);
            ws.dd(:, [ col.PP(1) col.PP(2)]) = ws.dd(:, [ col.PP(1) col.PP(2)])./repmat(sum(ws.dd(:, [ col.PP(1) col.PP(2)]),2),1,2);
        end
        
        subjdata{s, details.session+1}= ws.dd(:, [col.PL(1) col.PP(1) col.PL(2) col.PP(2) col.Choice col.ChoiceRT]);  % Order here must match fp.col
        details.subj_ntrials(s,1)=size(subjdata{s, details.session+1},1);
    end
    for m=startm:details.n_models 
        disp(['Model ' num2str(m) ' - ' details.models{m,1} '  ##############' ]); wm=[]; 
        r_grid{m,1}=details.models{m,1};
        
        % Set up parameters
        wm.startpar=cell(1,details.models{m,2});
        for p=1:details.models{m,2}
            x=details.models{m,6}{p};
            eval(['wm.startpar{p}=' details.models{m,5}{p} ';']);
        end
        details.models{m,7}=wm.startpar;
        
        for s=1: details.n_subjs  % Fit all subjects 
            disp(['Subject ' num2str(s) '  (' details.subjects{s} ')'])
            r_grid{m,2}{s,1}=details.subjects{s}; ws=[]; 
            ws.clockstart =clock;  
            
            % Execute 
            eval(['[ws.nLL, ws.np, ws.ns] = f_gridLCA(@' details.models{m,1} ', {subjdata{s, details.session+1} fp}, wm.startpar);'])
            r_grid{m,2}{s,3}=ws.nLL;
            ws.clockend=clock; 
            
            subfit=[]; 
            subfit.duration = [num2str(fmat(ws.clockend - ws.clockstart,4)) ' hr '  num2str(fmat(ws.clockend - ws.clockstart,5)) ' min ' num2str(fmat(ws.clockend - ws.clockstart,6),2) ' s']; 
            subfit.date= date;
            subfit.model= details.models{m,1}; 
            subfit.n_parsteps = details.n_parsteps;
            L=ws.nLL;  
            eval([  '[' strrep(strjoin(cellfun(@(x)[' ws.ind(' num2str(x) ')'],num2cell(1:details.models{m,2}), 'UniformOutput',0)),',,',',')  ']= ind2sub(size(L), find(L==min(L(:))));'])
            subfit.minL_indices=  ws.ind;  % Indices of best point 
             
            ws.filename= f_newname([details.subjects{s} '_se' num2str(details.session) 'stp' num2str(details.n_parsteps) 'sam' num2str(fp.nsamples)  '_mod' details.models{m,1} '_fit.mat'], [where.mod filesep '2 Inputs' filesep '3 Subgrids']);
            save([where.mod filesep '2 Inputs' filesep '3 Subgrids' filesep ws.filename], 'subfit', 'details', 'L') 
            disp( ['FIT   ' details.models{m,1} '  -  ' details.subjects{s}  ' [duration: '  subfit.duration '] -------------']) 
        end  
    end
    
    % Notify 
    try;  f_sendemail('kurzlich', ['[' details.whichPC '] Subject grids are done for model ' details.models{m,1} '  [' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ']'], ' ',1); end 
end

%% Assemble grid files into res files 
% Fit grid: 'r_grid' - model grids for each subject
%       Col 1=Model name
%       Col 2=Grid search results (nLL, parameters)
%             Col 1: Subject
%             Col 2: nLL
%             Col 3: Full grid of nLL (all param-points, p-dimensions)
%             Col 4: Indices of best fit (calculated later)
%             Col 6 onwards: best fit model parameters (1st parameter is beta/inverse temperature)

