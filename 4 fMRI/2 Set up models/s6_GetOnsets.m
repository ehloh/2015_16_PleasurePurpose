% Get onsets for set up of 1st level
% clear all; close all hidden; path(pathdef); clc;
clear all; path(pathdef); clc;

% Requested analysis 
logg.specificsubjects={
%     'p01_YH'; 'p02_MI'; 
%     'p03_AY';'p04_AW';'p05_CA';'p06_BT';'p07_HC';'p08_KC';'p09_KJ';'p10_YC';'p11_BI';'p12_AL';'p13_MS';'p14_SK';'p15_MM';'p16_SH';'p17_BB';'p18_WB';'p19_HB';'p20_LZ' 
    };

logg.func_prefix='s8wubf';
 
% Which onsets model? If scores are binned, this should say so in the model name ** ######################## 
% logg.onsetsmodel='r1b_PLnPPrating';  
logg.onsetsmodel='cl2b_ChoUnchoPLnPP';  
% logg.onsetsmodel='cl3b_PLnPP_conflict';  

for o1=1:1 % General settings and specifications 
    
    % Load subjects
    w=pwd; 
    if strcmp(w(1), '/')==1;  where.where='/Users/EleanorL/Dropbox/SCRIPPS/3b PLPP fmri';  where.data_brain='/Users/EleanorL/Desktop/3 PLPR/1 Brain data';   
        where.data_beh = '/Users/EleanorL/Dropbox/SCRIPPS/3 Pleasure purpose/3 Behaviour';
        where.behscripts = '/Users/EleanorL/Dropbox/SCRIPPS/3 Pleasure purpose/4a Beh analysis basic';
    else where.where='C:\Users\e.loh\Dropbox\SCRIPPS\3b PLPP fmri';  
        where.data_brain= 'D:\1 PLPP\1 MRI data';  where.spm='C:\toolbox\64\spm8'; 
        where.data_beh = 'C:\Users\e.loh\Dropbox\SCRIPPS\3 Pleasure purpose\3 Behaviour'; 
        where.behscripts = 'C:\Users\e.loh\Dropbox\SCRIPPS\3 Pleasure purpose\4a Beh analysis basic';
    end   
    [n t logg.datalog]=xlsread([where.data_brain fs  'datalog_plpr.xlsx']); 
    path(pathdef),  addpath(where.where),     addpath(where.behscripts),  addpath( [where.where filesep '2 Set up models' filesep '1 Onset models']);
    [logg.subjects logg.n_subjs logg.datalog] = f_selectsubjects(logg.datalog, logg.specificsubjects, [logg.datalog vertcat('include_all', num2cell(ones(size(logg.datalog,1)-1,1)))], 'include_all');
    if strcmp(logg.onsetsmodel(1:2), 'cl')  
        logg.specificsubjects= logg.subjects(~strcmp(logg.subjects, 'p12_AL'));  disp('Excluding p12_AL from this model');  
        [logg.subjects logg.n_subjs logg.datalog] = f_selectsubjects(logg.datalog, logg.specificsubjects, [logg.datalog vertcat('include_all', num2cell(ones(size(logg.datalog,1)-1,1)))], 'include_all');
    end 
    scan=load([where.where filesep '1 Preprocessing' fs 'i_scanningdetails.mat']);
    logg.FLthread=[]; f_mat=@(A,x)A(x);  errorlog=cell(1,1); e=1; 
    
    % Get columns 
	ws.f= load([where.data_beh fs logg.subjects{1} fs logg.subjects{1} '_behdata.mat']);
    col1= ws.f.col1; col2= ws.f.col2; col3= ws.f.col3; 

    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp(['No. of subjects: ' num2str(logg.n_subjs)])
    if isempty(logg.specificsubjects)==0; disp('   Subset of subjects only:'); disp(logg.subjects); end
    disp(['Data location: ' where.data_brain])
    disp('=======================================================')
    
end

 where.data_brain = '/Users/EleanorL/Dropbox/WorkPC/PLPP'; 

%% (1) Format data 

for o=1:1 % New columns  
    col1.PLbinrating = col1.PLbin;
    col1.PPbinrating = col1.PPbin;
    col1.PLPPbinrating = col1.PLPPbin;
end
 
disp('Fetching data + marking variables ############')
subjdata=[logg.subjects cell(logg.n_subjs,5)]; % 2=name, 3= rating session, 4= choice, 5= choice lab, all choice
d_sublog= cell(logg.n_subjs,1);  d_nscans = nan(logg.n_subjs,3);  % n scans in each
for s=1:logg.n_subjs % counting scans
    disp(['Subject ' num2str(s) '   -  ' logg.subjects{s}])
    
    % Load + sort behavioural data 
    ws.f= load([where.data_beh fs logg.subjects{s} fs logg.subjects{s} '_behdata.mat']);
    ws.d1 = ws.f.data1;   ws.d2 = ws.f.data2; ws.d3 = ws.f.data3; 
    ws.d2(:, col2.Session)=2;  ws.d3(:, col3.Session)=3; 
    
    % Check that preprocessing details are standard 
     d_sublog{s} = ws.f.log;  
     if s >1; 
         for i=1:length(d_sublog{s}.define.scorebins); if sum(d_sublog{s}.define.scorebins{i} - d_sublog{1}.define.scorebins{i})>0; error('Discrepancy in score bins!'); end ; end
         if strcmp( char(d_sublog{s}.fcf_raw), char(d_sublog{1}.fcf_raw))~=1 | strcmp( char(d_sublog{s}.fcf_binned), char(d_sublog{1}.fcf_binned))~=1 error('Discrepancy in conflict formulation!'); end   
     end 
         
     % Calculate new data here 
     ws.d1(:, col1.TrialOK)=  (~isnan( ws.d1(:, col1.PP))).*(~isnan( ws.d1(:, col1.PL)));
     ws.d2(:, col2.TrialOK)= floor( mean(~isnan( ws.d2(:, [col2.PL1 col2.PL2 col2.PP1 col2.PP2])),2)).*(~isnan(ws.d2(:, col2.Choice)));
     ws.d3(:, col3.TrialOK)= floor( mean(~isnan( ws.d3(:, [col3.PL1 col3.PL2 col3.PP1 col3.PP2])),2)).*(~isnan(ws.d3(:, col3.Choice)));
     ws.d2(:, col2.Duration_OptionsChoice)= ws.d2(:, col2.Onset_Choice)-ws.d2(:, col2.Onset_Options);
     ws.d2(:, col2.Duration_Trial)= ws.d2(:, col2.Onset_Choice)-ws.d2(:, col2.Onset_Options)+0.5;  % Assume end trial = 0.5s after choice. Fixation onsets are not marked in data
     ws.d2(:, col2.Duration_Options)=5; 
     ws.d3(:, col3.Duration_OptionsChoice)= ws.d3(:, col3.Onset_Choice)-ws.d3(:, col3.Onset_Options);
     ws.d3(:, col3.Duration_Trial)= ws.d3(:, col3.Onset_Choice)-ws.d3(:, col3.Onset_Options)+0.5;  % Assume end trial = 0.5s after choice. Fixation onsets are not marked in data
     ws.d3(:, col3.Duration_Options)=5; 
       
     % Output 
     subjdata{s,2}.d(1:3)= {ws.d1 ws.d2 ws.d3} ; 
    if strcmp(logg.onsetsmodel(1:3), 'ca_')  % If concat-choice requested, count scans + adjust onsets  
        if strcmp(logg.subjects{s}, 'p12_AL'), ws.func_prefix=strrep(logg.func_prefix, 'u','r'); else ws.func_prefix=logg.func_prefix; end
        ws.whereFL=[where.data_brain filesep logg.subjects{s} filesep '2 First level' filesep];
        d_nscans(s,1)=  size(spm_select('List', [ws.whereFL 'Preproc func s1'], ['^' ws.func_prefix '.*.img']) ,1);
        d_nscans(s,2)=  size(spm_select('List', [ws.whereFL 'Preproc func s2'], ['^' ws.func_prefix '.*.img']) ,1);
        d_nscans(s,3)=  size(spm_select('List', [ws.whereFL 'Preproc func s3'], ['^' ws.func_prefix '.*.img']) ,1);
        disp(['     No. Scans: ' num2str(d_nscans(s,:))])  
        ws.d4=[ws.d2; ws.d3];  ws.d4(ws.d4(:, col3.Session)==3, [col3.Onset_Options col3.Onset_Choice])= ws.d4(ws.d4(:, col3.Session)==3, [col3.Onset_Options col3.Onset_Choice]) + sum(d_nscans(s,2)) * (scan.TRms*scan.nSlicesPerVol/1000);
        subjdata{s,2}.d{4}=  ws.d4; 
    end 
    ws=[]; 
end
 
%% Checks and messing w data 
 

do_checks = 0;
if do_checks
    s=2;    disp(logg.subjects{s})
    d1= subjdata{s,2}.d{1};
    d2= subjdata{s,2}.d{2};
    d3= subjdata{s,2}.d{3};  
    d3=d3(d3(:, col3.TrialOK)==1, :);
    
    % After binning scores for p02 session 2, are options 1 & 2 exactly identical?
    dd=  d2(:,  [col2.choPL col2.choPP col2.unchoPL col2.unchoPP]); disp('looking at main choice session!')  
    dd=  d3(:,  [col3.choPL col3.choPP col3.unchoPL col3.unchoPP]); disp('looking at choice lab session!')  
    for t=1:size(dd,1)
        for tt=1:4
            if dd(t,tt)<3.5, ddb(t,tt)=1;
            elseif dd(t,tt)>7.5, ddb(t,tt)=3;
            else  ddb(t,tt)=2;
            end
        end 
    end
%     [ddb d3(:,  [col3.choPLbin col3.choPPbin col3.unchoPLbin col3.unchoPPbin])]
    
    % Sanity checks for Session 2, Chosen & unchosen logged properly?    
    if sum(nansum(abs(d2(d2(:, col2.Choice)==1, [col2.choPL col2.choPP col2.unchoPL col2.unchoPP])-d2(d2(:, col2.Choice)==1, [col2.PL1 col2.PP1 col2.PL2 col2.PP2 ])))) + sum(nansum(abs(d2(d2(:, col2.Choice)==2, [col2.choPL col2.choPP col2.unchoPL col2.unchoPP])- d2(d2(:, col2.Choice)==2, [col2.PL2 col2.PP2 col2.PL1 col2.PP1]))))>0; disp('Discrep in Cho vs Uncho assignments!'); end  
    for t=1:size(d2, 1)
        wt.ev1 = d2(t, [col2.Option1 col2.PL1 col2.PP1]); wt.ev2 = d2(t, [col2.Option2 col2.PL2 col2.PP2]); 
        % Session 2, Chosen & unchosen logged properly?
        if (sum(wt.ev1(2:3) -d1(d1(:,col1.Eventnum)==wt.ev1(1), [col1.PL col1.PP])) + sum(wt.ev2(2:3) -d1(d1(:,col1.Eventnum)==wt.ev2(1), [col1.PL col1.PP]))) >0;   disp( ['Discrep choice trial ' num2str(t)]); end
    end
end
%% (2) Format regressors
% % subjdata: Col 1= Name, Col 2=data {.d} (1Rating, 2Choice, 3ChoiceLab, 4ChoiceAll), Col3=onsets variables 
input('Format regressors?   ');

for i=1:logg.n_subjs
    s=find(strcmp(subjdata, logg.subjects{i})==1); c=1;
    disp(['Subject ' num2str(i) '  (' logg.subjects{s} ')'])
    
    % Create variables, according to requested model
    eval(['[ws.variables c] ='  logg.onsetsmodel   '(subjdata{s,2}.d, {col1 col2 col3}, c);']);
    
    subjdata{s,3}=ws.variables;
    ws=[];
end

%% (3) Save onsets ------------------------------
 


% Save  
disp('Saving onsets to First Level folders')
for i=1:logg.n_subjs
    s=find(strcmp(subjdata, logg.subjects{i})==1); % Identify correct subject
    names=subjdata{s,3}.names;
    onsets=subjdata{s,3}.onsets;
    durations=subjdata{s,3}.durations;
    if isfield(subjdata{s,3},  'pmod');  pmod=subjdata{s,3}.pmod;  else pmod=[];  end
    details =  d_sublog{s};
    if strcmp(f_mat(logg.onsetsmodel, f_mat(regexp(logg.onsetsmodel(1:f_mat(strfind(logg.onsetsmodel, '_'),1)), '\d'), length(regexp(logg.onsetsmodel(1:f_mat(strfind(logg.onsetsmodel, '_'),1)), '\d')))+1), 'b'); % Is this a bin modeled? Binned-score models have a b right after the number in the prefix **b_ (or bi_, bii_ .. i/ii referring to bin length),
        details.bin_nlevels  = length(details.define.scorebins);
    end
    
    for k=1:length( onsets),
         if (sum( isnan(onsets{k})) + sum( isinf(onsets{k})))>0 ; error(['WARNING: Bad onsets for ' logg.subjects{s} ' - ' names{k}]), end 
         if (sum( isnan(durations{k})) + sum( isinf(durations{k})))>0 ; error(['WARNING: Bad durations for ' logg.subjects{s} ' - ' names{k}]), end
         
         if ~isempty(pmod) &  k<=length(pmod) & ~isempty(pmod(k).name)
             for p=1:length(pmod(k)) 
                 if sum(isnan(pmod(k).param{p})) + sum(isinf(pmod(k).param{p}))~=0 
                     error(['WARNING: Bad pmod values for ' logg.subjects{s} ' -  ' pmod(k).name{p}  ' pmod on onsets ' names{k} ')'])  
                 end 
                 if length(unique(pmod(k).param{p} ))==1;  
                     error(['WARNING: Only a single pmod value for ' logg.subjects{s} ' -  ' pmod(k).name{p}  ' pmod on onsets ' names{k} ')']) 
                 end 
             end 
         end
    end     

    save([where.data_brain filesep subjdata{s,1} filesep '2 First level' logg.FLthread filesep subjdata{s,1} '_onsets_' logg.onsetsmodel '.mat'], 'names', 'onsets', 'durations', 'pmod','details');
end

%% END

disp('############################################################'), disp('END');  disp(' ')
disp('Errors:');  disp(errorlog)


%% CHECK ONSETS for nans etc
% Manually clear all and load onsets + execute following 


docheck=0;
if docheck
    input('Manually check onsets etc?'); 
    for k=1:length(durations)  % Onsets and durations
        if sum(isnan(durations{k})) + sum(isinf(durations{k}))==0 & sum(isnan(onsets{k})) + sum(isinf(onsets{k}))==0
            disp([names{k} ':  onsets and durations ok']);
        else  disp([names{k} ':  onsets and durations BAD ']);
        end
    end
    for k=1:length(pmod)  % pmod contents
        if isempty(pmod(k))==1
            for p=1:length(pmod(k))
                if sum(isnan(pmod(k).param{p})) + sum(isinf(pmod(k).param{p}))==0
                    disp([names{k} '    -   ' pmod(k).name '  :  pmod values ok']);
                else  disp([names{k} '    -   ' pmod(k).name '  :  pmod values BAD ']);
                end
            end
        elseif pmod(k).name
            onsets(k)
            pmod(k)
            onsets(:).name
        end
    end
end

