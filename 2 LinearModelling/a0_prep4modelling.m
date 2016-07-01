% Load data file (from basic beh analysis fol), append new data for modelling 
clear all;close all hidden; clc

% Request specific
logg.specificsubjects={}; % BLANK to process all subjects
request.split_ICCC=1;  

request.choicesession=2;    % 2= fMRI, 3= in lab 
request.dataset='All data b3 (17-Nov-2015)';
% request.dataset='All pilot data (11-Nov-2015)';request.choicesession=1 ;    
 
for o1=1:1 % General settings and specifications
   
    % Load subjects
    where.where='/Users/EleanorL/Dropbox/SCRIPPS/7 Pleasure purpose';   
    where.mod =[where.where filesep '4b Modelling'];  where.beh=[where.where filesep '3 Behaviour'];  addpath(where.where)
    if ischar(logg.specificsubjects)==1 && strcmp(logg.specificsubjects, 'AllDataOK'), logg.specificsubjects= {'p01_YH';'p02_MI';'p03_AY';'p04_AW';'p05_CA';'p06_BT';'p07_HC';'p08_KC';'p09_KJ';'p10_YC';'p11_BI';'p13_MS';'p14_SK';'p15_MM';'p16_SH';'p17_BB';'p18_WB';'p19_HB';'p20_LZ';}; end ;
    
    if isempty(strfind(request.dataset, 'pilot'))==1 & isempty(strfind(request.dataset, 'Pilot'))==1
        [n t r]=xlsread([where.beh filesep 'datalog_plpr.xlsx']);
        [logg.subjects logg.n_subjs logg.datalog] = f_selectsubjects(r, logg.specificsubjects,r, 'All');
    else
        logg.subjects=[cellfun(@(x)['p0' num2str(x)], num2cell(1:9),'UniformOutput',0)'; cellfun(@(x)['p' num2str(x)], num2cell(10:20),'UniformOutput',0)';];
        logg.n_subjs = length(logg.subjects); 
        [logg.subjects logg.n_subjs logg.datalog] = f_selectsubjects([[{'Subject'}; logg.subjects]  [{'All'}; num2cell(ones(20,1))] ], logg.specificsubjects,[[{'Subject'}; logg.subjects]  [{'All'}; num2cell(ones(20,1))] ], 'All');
    end
%     [n t r]=xlsread([where.beh filesep 'datalog_plpr.xlsx']); 
%     [logg.subjects logg.n_subjs logg.datalog] = f_selectsubjects(r, logg.specificsubjects,r, 'All'); 
    w= load( [where.mod fs '2 Inputs' fs  '2 Data' fs   request.dataset] );  subjdata=w.subjdata;  col=w.col;  
    try logg.define=w.log.define; catch;  try logg.define=w.logg.define;  catch; disp('Could not find log.define'); end;  end
    [d d subjdata] = f_selectsubjects([{'Subject'}  cell(1, size(subjdata,2)-1 )  ; subjdata ], logg.subjects,[[{'Subject'}  cell(1, size(subjdata,2)-1 ) {'All'} ]; [subjdata num2cell(ones(size(subjdata,1),1))]], 'All');  subjdata=subjdata(2:end, :);
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
    disp(['No. of subjects: ' num2str(logg.n_subjs)])
    if isempty(logg.specificsubjects)==0;  if logg.n_subjs == 19 & sum(strcmp(logg.subjects, {'p01_YH';'p02_MI';'p03_AY';'p04_AW';'p05_CA';'p06_BT';'p07_HC';'p08_KC';'p09_KJ';'p10_YC';'p11_BI';'p13_MS';'p14_SK';'p15_MM';'p16_SH';'p17_BB';'p18_WB';'p19_HB';'p20_LZ';}) ) ==19; disp('   Subjects w complete data only (n=19)'); else disp('   Subset of subjects only:'); disp(logg.specificsubjects); end;   end; disp(' ')
    disp(' '); input('Hit Enter to start      ')
    disp('=======================================================')
    
end

%% Define quantities 

for o=1:1 % Archive of un-used 
    f_rescale=@(rmin, rmax,x)    rmin + (rmax-rmin).*(x- min(x(:))) ./ max(x- min(x(:))); 
    
    % Note: conflict is already computed by this point. If you want to play w other formulations, 
    % give them other names. Note also that if you want to set up other IVs that are 
    % not choice invariant (i.e. dependent on chosen option's value, not just available options' values),
    % the RT analysis needs to be flipped.

end

%% Implementation

for o=1:1  % Data columns 
    
    % Session 2 
    w.s2.colmax=structmax(col.s2);
    
end
eval(['ccol=col.s' num2str(request.choicesession) ';'])  % Columns for computation, not saved 

% Implement changes 
for s=1:logg.n_subjs
    ws.d= subjdata{s,request.choicesession+1};
    %
    
    if request.split_ICCC
        ws.ic_trials=   sum([(ws.d(:, ccol.PLbest)==1 & ws.d(:, ccol.PPbest)==2)  (ws.d(:, ccol.PLbest)==2 & ws.d(:, ccol.PPbest)==1)],2);
        subjdata_ic{s,request.choicesession+1} = ws.d(find(ws.ic_trials),:); 
        subjdata_cc{s,request.choicesession+1} = ws.d(find(1-ws.ic_trials),:); 
    end
    
    %
    subjdata{s,request.choicesession+1} = ws.d;
    ws=[]; 
end


% Note: Sessions 1 & 3 are not set up at all yet.

%% End

logg.col=col; logg.dataset= request.dataset;
save( [where.mod fs '2 Inputs' fs '2 Data' fs    request.dataset], 'subjdata', 'logg','col' );
if request.split_ICCC
    subjdata=subjdata_ic;  subjdata(:, 1)=logg.subjects;
    save( [where.mod fs '2 Inputs' fs '2 Data' fs  request.dataset(1:strfind(request.dataset, '(')-1) 'IC '  request.dataset(strfind(request.dataset, '('):end) ], 'subjdata', 'logg','col' );
    subjdata=subjdata_cc; subjdata(:, 1)=logg.subjects;
    save( [where.mod fs '2 Inputs' fs '2 Data' fs  request.dataset(1:strfind(request.dataset, '(')-1) 'CC '  request.dataset(strfind(request.dataset, '('):end) ], 'subjdata', 'logg','col' );
end 
    

if isempty(logg.define)==0 
disp('########################'), disp('[Quantifications:]');  
disp(' '),  

disp([char(logg.define(:,1)) char(repmat({'  '}, size(logg.define,1), 1)) char(logg.define(:,2))] )
end
