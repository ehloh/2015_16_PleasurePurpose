% First level analysis - Specify & Estimate model
clear all; path(pathdef); clc;

% Requested analysis 
log.specificsubjects={
%     'p01_YH';
'p02_MI';
%     'p03_AY';'p04_AW';'p05_CA';'p06_BT';'p07_HC';'p08_KC';'p09_KJ';'p10_YC';'p11_BI';'p12_AL';'p13_MS';'p14_SK';'p15_MM';'p16_SH';'p17_BB';'p18_WB';'p19_HB';'p20_LZ' 
    };

% Request steps 
request.specify=1;
request.RemoveEventOnsetsBeforeEstim=1; % Value, Competing & Orthog onsets families only
request.estimate=1;
 
% Which onsets model? ########################      
% log.onsetsmodel='r1b_PLnPPrating';  
log.onsetsmodel='cl2b_ChoUnchoPLnPP';  
% log.onsetsmodel='cl3b_PLnPP_conflict';  

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
    [n t log.datalog]=xlsread([where.data_brain fs  'datalog_plpr.xlsx']); 
    path(pathdef),  addpath(where.where),     addpath(where.behscripts),  addpath( [where.where filesep '2 Set up models' filesep '1 Onset models']);
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    if strcmp(log.onsetsmodel(1:2), 'cl')  
        log.specificsubjects= log.subjects(~strcmp(log.subjects, 'p12_AL'));  disp('Excluding p12_AL from this model');  
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    end 
    
    % Stable settings 
    request.includederivatives=1; log.FirstLevelThread=[]; log.prefix='s8wubf'; 
    scan=load([where.where filesep '1 Preprocessing' fs 'i_scanningdetails.mat']);
    f_mat=@(A,x)A(x); errorlog=cell(0,1); e=1; 
    
    % What sort of model is this?  
    switch log.onsetsmodel(1: regexp(log.onsetsmodel, '\d')-1)
        case 'r',  log.modeltype= 'Rating';  
        case 'c',  log.modeltype= 'Choice';  
        case 'cl',  log.modeltype= 'ChoiceLab';  
        otherwise  error('Scripts are not yet set up for this type of fMRI model yet'); 
    end  
    

    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
    disp(['Data location: ' where.data_brain])
    disp(['FL model: '  log.onsetsmodel])
    disp('=======================================================')
    
end
 
%% STEP 1:  Specify model
 
if request.specify==1; % Specify: Format trial onsets & other regressors
    disp('STEP 1: Specifying model ################################')
    for o1=1:1 % Settings 
        settings.firstlevelmodelspec.timing.units = 'secs';
        settings.firstlevelmodelspec.timing.RT = scan.TRms/1000*scan.nSlicesPerVol;
        settings.firstlevelmodelspec.timing.fmri_t = 16;
        settings.firstlevelmodelspec.timing.fmri_t0 = 1;
        settings.firstlevelmodelspec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
        settings.firstlevelmodelspec.sess.regress = struct('name', {}, 'val', {});
        settings.firstlevelmodelspec.sess.multi_reg = {''};
        settings.firstlevelmodelspec.sess.hpf = 128;
        settings.firstlevelmodelspec.fact = struct('name', {}, 'levels', {});
        if request.includederivatives==1; settings.firstlevelmodelspec.bases.hrf.derivs =[1 1];
        else settings.firstlevelmodelspec.bases.hrf.derivs =[0 0]; disp('Derivatives are NOT included');
        end
        settings.firstlevelmodelspec.volt = 1;
        settings.firstlevelmodelspec.global = 'None';
        settings.firstlevelmodelspec.mask = {''};
        settings.firstlevelmodelspec.cvi = 'AR(1)';
%         settings.firstlevelmodelspec.cond.sess.pmod.poly = 1; % Polynomials?
    end
    for s=1: log.n_subjs % Specify model for each subject
        disp(['Subject ' num2str(s) '   (' log.subjects{s} ') -------------------'])
        matlabbatch{1}.spm.stats.fmri_spec= settings.firstlevelmodelspec;
        try
            wb.where=[where.data_brain filesep log.subjects{s} filesep '2 First level' log.FirstLevelThread filesep];
             wb.wheremodel=[wb.where log.onsetsmodel ' Estimated' filesep]; if isdir(wb.wheremodel)==0; mkdir(wb.wheremodel); end
            wb.prefix = log.prefix;  if strcmp(log.subjects{s}, 'p12_AL');  wb.prefix =strrep(wb.prefix , 'u', 'r'); end 
             
            % SPECIFY MODEL -------------------  
            matlabbatch{1}.spm.stats.fmri_spec.dir = {wb.wheremodel}; % Specification files
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi ={[wb.where log.subjects{s} '_onsets_' log.onsetsmodel '.mat']};
            if exist([wb.where log.subjects{s} '_onsets_' log.onsetsmodel '.mat'], 'file')==0; input(['ERROR: Could not find onsets for the specified model   --   ' log.onsetsmodel ]); end
            
            % 
            switch log.modeltype
                case 'Rating';
                    wb.wherefuncs = [wb.where 'Preproc func s1' filesep ];  
                    f=spm_select('List', wb.wherefuncs , ['^' wb.prefix '.*img$']);  % This does assume functional scans are in order of runs!
                    wb.func=cellstr([repmat([wb.wherefuncs], size(f,1),1) f]);
                    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {[wb.where log.subjects{s} '_1rating_regmovement.txt']};
                case 'Choice' 
                    wb.wherefuncs = [wb.where 'Preproc func s2' filesep ];  
                    f=spm_select('List', wb.wherefuncs , ['^' wb.prefix '.*img$']);  % This does assume functional scans are in order of runs!
                    wb.func=cellstr([repmat([wb.wherefuncs], size(f,1),1) f]);
                    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {[wb.where log.subjects{s} '_2choice_regmovement.txt']}; 
                case 'ChoiceLab'
                    wb.wherefuncs = [wb.where 'Preproc func s3' filesep ];
                    f=spm_select('List', wb.wherefuncs , ['^' wb.prefix '.*img$']);  % This does assume functional scans are in order of runs!
                    wb.func=cellstr([repmat([wb.wherefuncs], size(f,1),1) f]);
                    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {[wb.where log.subjects{s} '_3choicelab_regmovement.txt']};
                otherwise, error('Not set up yet for this model type') %                     p02_MI_2choice_regmovement %                     p02_MI_3choicelab_regmovement   %                     p02_MI_4choiceall_regmovement
            end 
            
            
            % Choose functional files
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans =cellstr(wb.func);     %   save('batch','matlabbatch')
            spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
            if strcmp('b', log.onsetsmodel(regexp(log.onsetsmodel(1:f_mat(strfind(log.onsetsmodel, '_'),1)), '\d')+1))  % Have scores been binned for this model?
                wb.so = load([wb.where log.subjects{s} '_onsets_' log.onsetsmodel '.mat']);
                wb.s =load([wb.wheremodel 'SPM.mat']); SPM = wb.s.SPM;
                SPM.details.scorebins =  {wb.so.details.bin_nlevels wb.so.details.define.scorebins};  % Save scorebin details to SPM variable
            end
            matlabbatch=[];wb=[];  SPM=[];
        catch;  errorlog{e}=['ERROR: Could not specify model for subject   ' log.subjects{s} ]; disp(errorlog{e}); e=e+1;
        end
        ws=[];
    end
end

%% STEP 2: Remove unecessary regressors

if request.RemoveEventOnsetsBeforeEstim==1;
    disp('STEP 2: Removing regressors before model estimation ##############')
 
    % Auto-identify pmods (assumes all of-interest regresssors are identical for all subjects)
    %       Assume event onsets (regardless of session) are prefixed 'o_'
    request.regremove_prefix='o_';
    ws.where=[where.data_brain filesep log.subjects{1} filesep '2 First level' log.FirstLevelThread filesep log.onsetsmodel ' Estimated' filesep];
    ws=load([ws.where 'SPM.mat']);    ws.regnames=ws.SPM.xX.name(1:3:end)';   % Identifying 1st deriv only 
    request.regs2remove =      f_mat( ws.regnames(cellfun(@(x)~isempty(x), strfind(ws.regnames,  ['Sn(1) ' request.regremove_prefix]))) ,    find(cellfun(@(x)isempty(x), strfind(ws.regnames(cellfun(@(x)~isempty(x), strfind(ws.regnames,  ['Sn(1) ' request.regremove_prefix]))), '^'))));
   
    disp(['Regressors to remove (spare onsets assumed prefixed= ' request.regremove_prefix  '):'])
    disp(request.regs2remove), disp('--------------------------------------------------------------------------------')
     
    
    % Remove regressors
    for s=1:log.n_subjs         
        disp(['Subject ' num2str(s) '   (' log.subjects{s} ') -------------------'])
        try
            % Load details from the model to be altered
            wb.where=[where.data_brain filesep log.subjects{s} filesep '2 First level' log.FirstLevelThread filesep log.onsetsmodel ' Estimated' filesep];
            ws=load([wb.where 'SPM.mat']);
            ws.regnames=ws.SPM.xX.name';
            
            % Compile list of regressors
            [ ws.regdetails ] = f_GetRegressorNums(ws.regnames,request.regs2remove);
            
            j=1;
            for i=1:size(ws.regdetails,1) % Expand list to include derivatives
                for r=1:length(ws.regdetails{i,3})
                    ws.removelist{j,1}=ws.regdetails{i,3}{r}; j=j+1;
                    ws.removelist{j,1}=[ws.regdetails{i,3}{r}(1:length(ws.regdetails{i,3}{r})-2) '2)']; j=j+1; % 2nd derivative
                    ws.removelist{j,1}=[ws.regdetails{i,3}{r}(1:length(ws.regdetails{i,3}{r})-2) '3)']; j=j+1; % 3rd derivative
                end
            end
            for i=1:length(ws.removelist)  % Check: Any requested removals that don't exist?
                if isempty(strfind(ws.regnames,ws.removelist{i}))
                    error('Error: Regressor scheduled for removal does not exist')
                end
            end
            
            % Remove regressors
            ws.newSPM = f_RemoveRegressor(ws.SPM,ws.removelist);
            
            % Save
            SPM=ws.newSPM;
            SPM.originalSPM=ws.SPM;
            SPM.note='SPM file modified before model specification.';
            save([wb.where 'SPM.mat'], 'SPM');
            ws=[]; SPM=[];
            
        catch;  errorlog{e}=['ERROR: Could not remove regs b4 model est for subject   ' log.subjects{s} ]; disp(errorlog{e}); e=e+1;
        end
    end
    
end

%% STEP 3: Estimate model

if request.estimate==1
    disp('STEP 3: Estimate model ################################')
    for s=1:log.n_subjs
        try
            disp(['Subject ' num2str(s) '   (' log.subjects{s} ') -------------------'])
            wb.where=[where.data_brain filesep log.subjects{s} filesep '2 First level' log.FirstLevelThread filesep log.onsetsmodel ' Estimated'];
            %
            f   = spm_select('List', wb.where, 'SPM.mat');
            matlabbatch{1}.spm.stats.fmri_est.spmmat = cellstr([wb.where filesep f]);
            matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
            %
            spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
            matlabbatch=[];wb=[];
        catch
            errorlog{e}=['ERROR: Could not estimate model for subject   ' log.subjects{s} ]; disp(errorlog{e}); e=e+1;
            wb=[];
        end
    end
end

%% END5

disp('======================================================='); w.c=clock; disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
disp('Analysis completeC:\Users\e.loh\'); disp(request)
disp(['No. of subjects: ' num2str(log.n_subjs)]); disp(' '); 
disp(log); disp('Errors:'); disp(errorlog'); disp(' ')
disp('=======================================================')
 
diary off
try % Notify researcher
    f_sendemail('kurzlich', strcat(['Analysis batchscript is complete (s7_Firstlevel_1SpecifyEstimate  -  '  log.onsetsmodel  ')']), ' ',1);
end
