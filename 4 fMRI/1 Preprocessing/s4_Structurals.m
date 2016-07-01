% Process structurals and misc 
clear all;close all hidden; clc
 

STOPPED HERE. Need to redo the mean structural, one subject is outof wack.


% Requested analysis 
request.meanstruct=1; 
log.specificsubjects={}; % BLANK to process all subjects

log.specificsubjects={     
%     'p01_YH';  'p02_MI'; 
% 'p03_AY';
% 'p04_AW'; 'p05_CA';'p06_BT';'p07_HC';'p08_KC'; 'p09_KJ'; 'p10_YC'; 
%     'p11_BI'; 
% 'p12_AL'; 
% 'p13_MS';'p14_SK';'p15_MM';'p16_SH';'p17_BB';'p18_WB'; 'p19_HB';'p20_LZ'
}; % BLANK to request all subjects



for o1=1:1 % General settings and specifications
    
    % Load subjects
    w=pwd;  if strcmp(w(1), '/')==1;  where.where='/Users/EleanorL/Dropbox/SCRIPPS/3b PLPP fmri';  where.data_brain='/Users/EleanorL/Desktop/3 PLPR/1 Brain data';  
    
    where.spm='/Users/EleanorL/Documents/MATLAB/spm8';
    else where.where='C:\Users\e.loh\Dropbox\SCRIPPS\3b PLPP fmri';  where.data_brain= 'D:\1 PLPP\1 MRI data'; 
        where.spm='C:\toolbox\64\spm8';
    end   
     
    [n t log.datalog]=xlsread([where.data_brain fs  'datalog_plpr.xlsx']); addpath(where.where)
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    col.funcsess=[7 8 9];  col.funcfieldmaps=[10 11 12];  col.struct=13;  col.func_nscans=[14 15 16];  % log.datalog contents 
      
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp('Requested analysis:'); disp(request)
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.specificsubjects); end
    disp(' '); disp(['Data location: ' where.data_brain]); disp(' ')
    input('Hit Enter to start      ')
    disp('=======================================================')
    spm_jobman('initcfg');  f_subindex = @(A,r)A(r); 
end

%% (3) Calculate average structurals (T1w & MTw)

if request.meanstruct==1
    disp('############ Constructing average structurals ###############') 
     
    % Normalize individual subject structurals
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.template = { [where.spm filesep 'templates/T1.nii,1']};
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.weight = '';
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.smosrc = 8;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.smoref = 0;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.regtype = 'mni';
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.cutoff = 25;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.nits = 16;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = 1;
    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.preserve = 0;
    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.bb = [-78 -112 -50; 78 76 85];
    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.vox = [2 2 2];
    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.interp = 1;
    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.prefix = 'w';  
    for s=1:log.n_subjs
        ws.where=[where.data_brain filesep log.subjects{s} filesep '1 Preprocessed' filesep 'Structural' filesep];
        f=spm_select('List', ws.where, '^s.*img$'); 
        if isempty(f); error(['Cannot find structural for ' log.subjects{s} ]); end 
        
        matlabbatch{1}.spm.spatial.normalise.estwrite.subj(s).wtsrc = '';
        matlabbatch{1}.spm.spatial.normalise.estwrite.subj(s).resample =  {[ws.where f ',1']}; 
        matlabbatch{1}.spm.spatial.normalise.estwrite.subj(s).source =  {[ws.where f ',1']};  
    end 
    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    matlabbatch=[]; 
    
    % Constrcut group average
    disp(['[Group] IMCALC-ING GROUP SCAN----------'])
    matlabbatch{1}.spm.util.imcalc.output =  ['MeanStruc_n' num2str(log.n_subjs) '.img'];
    matlabbatch{1}.spm.util.imcalc.outdir = {where.data_brain};
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
    matlabbatch{1}.spm.util.imcalc.expression='(';
    for s=1:length(log.subjects) % Generate expression
        matlabbatch{1}.spm.util.imcalc.expression=[matlabbatch{1}.spm.util.imcalc.expression 'i' num2str(s)];
        if s<log.n_subjs matlabbatch{1}.spm.util.imcalc.expression=[matlabbatch{1}.spm.util.imcalc.expression '+'];
        else matlabbatch{1}.spm.util.imcalc.expression=[matlabbatch{1}.spm.util.imcalc.expression ')/' num2str(log.n_subjs)];
        end
    end
    for s=1:length(log.subjects) % Collect scans
        ws.where=[where.data_brain filesep log.subjects{s} filesep '1 Preprocessed' filesep 'Structural' filesep];
        f=spm_select('List', ws.where, '^ws.*img$');   
        wb.t1s{s} =[ws.where f ',1'];
    end
    matlabbatch{1}.spm.util.imcalc.input =wb.t1s;
    %
    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    matlabbatch=[];
end

%% END

disp('====================================')
w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
disp('Analysis completed:')
disp(request); disp(['No. of subjects: ' num2str(log.n_subjs)]); disp(' ')
disp('====================================')

diary off
try % Notify researcher
    f_sendemail('kurzlich', strcat('Analysis batchscript is complete (s2_MPM)'), ' ',1);
catch
end
