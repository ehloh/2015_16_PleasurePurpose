% Second level analysis - Parametized (Compete/Orthog/Trialtype) models only
clear all;close all hidden; clc

% Requested analysis
log.specificsubjects={
%     'p01_YH';'p02_MI'; 'p03_AY';'p04_AW';'p05_CA';'p06_BT';'p07_HC';'p08_KC';'p09_KJ';'p10_YC';'p11_BI';'p12_AL';'p13_MS';'p14_SK';'p15_MM';'p16_SH';'p17_BB';'p18_WB';'p19_HB';'p20_LZ' 
    };
 

% Which model? ########################
% log.onsetsmodel='r1_PLnPPrating'; 
log.onsetsmodel='c2_ChoUnchoPLnPP';  
% log.onsetsmodel='cl3b_PLnPP_conflict';  




% Subsample ########################
log.subsample ={}; 
% log.subsample = {'regr_oe_PLvPP' 'Low'}; 
% log.subsample = {'regr_oe_PLvPP_abs_r04' 'All'};   
log.subsample = {'regr_ot_PLvPP_abs_r04' 'All'};   
% log.subsample = {'beh1.PLvPP_abs_r04' 'All'};


log.secondlevelmodel='CondVar_2x2';                 % SECOND LEVEL MODEL###############
% log.secondlevelmodel='Identity_1samplettest';
% log.secondlevelmodel='Identity_pairedttest'; 

for o1=1:1 % General settings and specifications   
    
    % Load subjects
    w=pwd; 
    if strcmp(w(1), '/')==1;  where.where='/Users/EleanorL/Dropbox/SCRIPPS/3b PLPP fmri';  
        where.experiment_folder = '/Users/EleanorL/Desktop/3 PLPR';
        where.data_brain='/Users/EleanorL/Desktop/3 PLPR/1 Brain data';   
        where.data_beh = '/Users/EleanorL/Dropbox/SCRIPPS/3 Pleasure purpose/3 Behaviour';
        where.behscripts = '/Users/EleanorL/Dropbox/SCRIPPS/3 Pleasure purpose/4a Beh analysis basic';
    else where.where='C:\Users\e.loh\Dropbox\SCRIPPS\3b PLPP fmri';  
        where.experiment_folder ='D:\1 PLPP'; 
        where.data_brain= 'D:\1 PLPP\1 MRI data';  where.spm='C:\toolbox\64\spm8'; 
        where.data_beh = 'C:\Users\e.loh\Dropbox\SCRIPPS\3 Pleasure purpose\3 Behaviour'; 
        where.behscripts = 'C:\Users\e.loh\Dropbox\SCRIPPS\3 Pleasure purpose\4a Beh analysis basic';
    end   
    [n t log.datalog]=xlsread([where.data_brain fs  'datalog_plpr.xlsx']); 
    path(pathdef),  addpath(where.where),     addpath(where.behscripts),  addpath( [where.where filesep '2 Set up models' filesep '2 SL models']);
    
    % Stable settings  
    log.FirstLevelThread=[]; log.prefix='s8wubf'; 
    scan=load([where.where filesep '1 Preprocessing' fs 'i_scanningdetails.mat']);
    f_mat=@(A,x)A(x); errorlog=cell(0,1); e=1; 
    
    % What sort of model is this? 
    switch log.onsetsmodel(1: regexp(log.onsetsmodel, '\d')-1)
        case 'r',  log.modeltype= 'Rating';  
        case 'c',  log.modeltype= 'Choice';  
        case 'cl',  log.modeltype= 'ChoiceLab';  
        otherwise,  error('Scripts are not yet set up for this type of fMRI model yet'); 
    end   
 
    
    % Subsample? 
    if isempty(log.specificsubjects) &&  ~isempty(log.subsample)
        log.specificsubjects= f_subsample(log.onsetsmodel, log.subsample{1});  
        
        switch log.subsample{2}
            case 'Low';  log.specificsubjects = log.specificsubjects(1:round(length(log.specificsubjects)/2));
            case 'All';
            otherwise, error('Sorting for this subsample now done yet!')
        end 
          
        % Name this subset thread (ignore dots in subsample)
        log.subsample_name = ['Sample ' strrep(log.subsample{1},'.', ' ') ' ' log.subsample{2}]; 
    elseif ~isempty(log.specificsubjects) && ~isempty(log.subsample)
        input('Subsample turned off because specific subjects requested. OK?  ');
    else log.subsample_name=[]; 
    end  
    
    % Subjects 
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    if strcmp(log.onsetsmodel(1:2), 'cl')
        log.specificsubjects= log.subjects(~strcmp(log.subjects, 'p12_AL'));  disp('Excluding p12_AL from this model');  
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    end  
    
    

    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0;disp(['   Subset of subjects only (n=' num2str(log.n_subjs) '):']); disp(log.subjects);  end
    disp(' '); disp(['Data location (brain): ' where.data_brain])
    disp(' ');disp('CHECK HERE: 1st and 2nd level models ---------'); disp(' ')
    disp(['             First level model:  ' log.onsetsmodel])
    disp(['             Second level model:  ' log.secondlevelmodel]); disp(' ')
    input('Hit Enter to start      ')
    disp('=======================================================')
    
end

%% Which Variables are included in this model? (sample 1st subject)

for o=1:1 % Parse FL model details 
    
    % Read out conditions of interest in this fMRI model (Assuming all subject have the same of-interest regressors)
    wc=load([where.data_brain filesep log.subjects{1} filesep '2 First level' log.FirstLevelThread filesep log.onsetsmodel ' Contrasted' filesep 'SPM.mat']);
    log.FLregs =  f_mat(cellfun(@(x)x(7:length(x)), wc.SPM.xX.name(1:3:end)','uniformoutput',0), 1:find(cellfun(@(x)strcmp(x(1:2), 'n_'), cellfun(@(x)x(7:length(x)), wc.SPM.xX.name(1:3:end)','uniformoutput',0)),1,'first')-1);   % Up to b4 1st no-interest regressor    
    log.FLcons =  cellstr(char(wc.SPM.xCon(:).name));   if isempty(log.FLcons ), error('No FL cons for this model (looked at subject 1)'); end 
     
    % Identify variables that have PL & PP as different conditions (replaced by marker @@@)
    wc.plppcons=   f_mat(log.FLcons , [find(cell2mat(cellfun(@(x)~isempty(strfind(x, 'PL')), log.FLcons,'uniformoutput',0)));  find(cell2mat(cellfun(@(x)~isempty(strfind(x, 'PP')), log.FLcons,'uniformoutput',0)))]); 
    log.PLPPvars =  cellfun(@(x)strrep(x, 'PP','@@@'), f_mat(strrep(wc.plppcons(find(cellfun(@(x)~isempty(strfind(x, 'PL')), wc.plppcons))), 'PL','PP'), find(cell2mat(cellfun(@(x,b)sum(strcmp(x,b)), strrep(wc.plppcons(find(cellfun(@(x)~isempty(strfind(x, 'PL')), wc.plppcons))), 'PL','PP'),  repmat({wc.plppcons}, length(strrep(wc.plppcons(find(cellfun(@(x)~isempty(strfind(x, 'PL')), wc.plppcons))), 'PL','PP')),1), 'uniformoutput',0) ))),'uniformoutput',0) ;  % Find variables that have both PL & PP components. One line of code bitch.
    
    disp('#################################################')
    disp('FL cons in this model:'); disp(' ');disp(log.FLcons) 
    disp('PL vs PP conditions/variables in this model:')
    disp(' '); disp(log.PLPPvars ); disp(' ');
    disp('#################################################')
  
    
end

% Request specific (empty to auto parse) ------------------------
%       CondVar_2x2:    See script documentation 
%       Identity_pairedttest:  {model/contrast name, {Cond1, Cond2}}   e.g. {'oe_@@@rating', {'PL', 'PP'}: @@@ marker is replaced by PL & PP to search for the contrasts to compare w each other 
%       Identity_1samplettest:  {model/contrast name}   e.g. {'oe_event', 'trial'} ... 1sample ttest will be performed for all requested
%   You can also keyboard the SL models themselves and change settings there! 

req.SL ={}; 
% req.SL ={'trial'}; 

for o=1:1   % Settings for 2nd level models 
    if isempty(req.SL)
        % Auto-parse regressors and conditions
        if  strcmp(log.secondlevelmodel, 'Identity_pairedttest')
            req.SL = {log.PLPPvars {'PL', 'PP'}};  % Search term, replacement terms (conditions)
            if length(log.PLPPvars)>1;  error('Paired ttest set up for 1 case. havent yet extended to the 2 varialbe scenariuo!'); end
        elseif  strcmp(log.secondlevelmodel, 'CondVar_2x2')
            % Assume it's PL/PP x Variable
            if length(unique(cellfun( @(x) x(1:strfind(x,'_')),  log.PLPPvars, 'uniformoutput',0)) )==1   % If there's only 1 time epoch, assume the variable relates to this (i.e. time epoch is not the 2x2 variable)
                req.SL{1} = ['PLvPP ' strrep( strjoin( cellfun(@(x)[upper(x(1)) x(2:end)], cellfun( @(x) strrep(x,unique(cellfun( @(x) x(1:strfind(x,'_')),  log.PLPPvars, 'uniformoutput',0)),''),  cellfun( @(x) strrep(x,'@@@',''),  log.PLPPvars, 'uniformoutput',0)), 'uniformoutput',0)' ), ',', '')];
                req.SL{2} =  [{'PL' 'PP'}; cellfun( @(x) strrep(x,unique(cellfun( @(x) x(1:strfind(x,'_')),  log.PLPPvars, 'uniformoutput',0)),''),  cellfun( @(x) strrep(x,'@@@',''),  log.PLPPvars, 'uniformoutput',0))'];
                req.SL{3} =  [cellfun( @(x) strrep(x,'@@@','PL'),  log.PLPPvars , 'uniformoutput',0)   {[1 1]; [1 2]}; cellfun( @(x) strrep(x,'@@@','PP'),  log.PLPPvars , 'uniformoutput',0)   {[2 1]; [2 2]}];
                disp([req.SL{1} ' SL model (PL vs PP is assumed to be one of the factors) --------------------------------- '])
                disp(req.SL{2}), disp(req.SL{3})
            else error('Requesting SL not set up yet')
            end
        else error('Requesting SL not set up yet')
        end
    end
end
 

%% (2) Second-level model specification + Estimation (In function)

disp('STEP 1: Specify 2nd-level model  ####################')

% Results thread for this first level model
where.resultsfolder=[where.experiment_folder filesep '2 Second level results' log.FirstLevelThread filesep log.onsetsmodel filesep log.subsample_name];   if isdir(where.resultsfolder)==0; mkdir(where.resultsfolder); end

input('USER: Continue to specify + estim model?     ');

% Specify requested model
eval([' [ matlabbatch contrastfiles] = ' log.secondlevelmodel '(log, where,  log.subjects, log.onsetsmodel, log.secondlevelmodel, req.SL);'])


%% END

disp('======================================================='); w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
disp('Analysis completed (s8_Secondlevel)'); disp(['No. of subjects: ' num2str(log.n_subjs)])
if isempty(log.specificsubjects)==0; disp('   Subset of log.specificsubjects only:'); disp(log.subjects); end
disp(' '); disp(['Data location (brain): ' where.data_brain]); disp(' ')
disp(['First level model:  ' log.onsetsmodel])
disp(['Second level model:  ' log.secondlevelmodel])
disp(' '); disp('Errors:'); disp(errorlog'); disp(' ')
disp('=======================================================')

diary off
try % Notify researcher
    f_sendemail('kurzlich', strcat('Analysis batchscript is complete (s8_Secondlevel)'), ' ',1);
end
cd(where.resultsfolder)