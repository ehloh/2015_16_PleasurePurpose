% First level contrasts - Parameterized (Compete/Orthog) & Trial type onsets models
clear all; path(pathdef); clc;

% Requested analysis 
log.specificsubjects={
%     'p01_YH'
'p02_MI'; 
% 'p03_AY';'p04_AW';'p05_CA';
%     'p06_BT';'p07_HC';'p08_KC';'p09_KJ';'p10_YC';'p11_BI';'p12_AL';'p13_MS';'p14_SK';'p15_MM';'p16_SH';'p17_BB';'p18_WB';'p19_HB';'p20_LZ' 
    };

% Request steps 
req.NonIdentityCons = 0; 
req.DeleteExistingFLContrasts=1;

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
    req.includederivatives=1; log.FirstLevelThread=[]; log.prefix='s8wubf'; 
    scan=load([where.where filesep '1 Preprocessing' fs 'i_scanningdetails.mat']);
    f_mat=@(A,x)A(x); errorlog=cell(0,1); e=1; 
    
    % What sort of model is this? 
    switch log.onsetsmodel(1: regexp(log.onsetsmodel, '\d')-1)
        case 'r',  log.modeltype= 'Rating';  
        case 'c',  log.modeltype= 'Choice';
        case 'cl',  log.modeltype= 'ChoiceLab';  
        otherwise;  error('Scripts are not yet set up for this type of fMRI model yet'); 
    end  
    

    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
    disp(['Data location: ' where.data_brain])
    disp('=======================================================')
    
end

%% Set up Contrasts
%     req.cons: {1} contrast name {2} {regname, regtype weight; ..} ;
%                    (regtype: 1=event, 2=pmod)

req.cons= cell(0,2); 

% Identity contrasts (all of-interest regressors, auto-identified ) 
%   Assume execute contrasts for all of-interest (~ prefix n_) regressors. Assumed that all of-interest regressors are identical for all subjects 
try wc=load([where.data_brain filesep log.subjects{1} filesep '2 First level' log.FirstLevelThread filesep log.onsetsmodel ' Estimated' filesep 'SPM.mat']);    
catch;  wc=load([where.data_brain filesep log.subjects{1} filesep '2 First level' log.FirstLevelThread filesep log.onsetsmodel ' Contrasted' filesep 'SPM.mat']);    
end; log.regnames=wc.SPM.xX.name(1:3:end)';   % Identifying 1st deriv only
req.iconregs =   cellfun(@(x)x(7: strfind(x, '*bf(1)')-1), log.regnames( 1: find(cellfun(@(x)~isempty(x), strfind(log.regnames, 'Sn(1) n_')),1)-1),'uniformoutput',0);   % Weight all regs up to (no including) first n_ prefix regrerssors
for i= 1:length (req.iconregs ) % Auto-read names of 
    if isempty(strfind(req.iconregs{i}, '^1')), req.cons(size(req.cons,1)+1, 1:2)  = {req.iconregs{i}    {req.iconregs{i} 1 1}};  
    elseif strcmp(log.modeltype,'Rating') 
        wc.conname= strrep( req.iconregs{i}(1:length(req.iconregs{i})-2), 'o_eventx', 'oe_');
        req.cons(size(req.cons,1)+1, 1:2)  = {wc.conname   {req.iconregs{i}(1:length(req.iconregs{i})-2) 2 1}};  
    elseif strcmp(log.modeltype,'Choice') 
        wc.conname= strrep( req.iconregs{i}(1:length(req.iconregs{i})-2), 'o_trialx', 'ot_');
        wc.conname= strrep( wc.conname, 'o_optionsx', 'op_');
        wc.conname= strrep( wc.conname, 'o_optionschoicex', 'oq_');
        req.cons(size(req.cons,1)+1, 1:2)  = {wc.conname   {req.iconregs{i}(1:length(req.iconregs{i})-2) 2 1}};   
    elseif strcmp(log.modeltype,'ChoiceLab')  
        wc.conname= strrep( req.iconregs{i}(1:length(req.iconregs{i})-2), 'o_trialx', 'ot_');
        wc.conname= strrep( wc.conname, 'o_optionsx', 'op_');
        wc.conname= strrep( wc.conname, 'o_optionschoicex', 'oq_');
        req.cons(size(req.cons,1)+1, 1:2)  = {wc.conname   {req.iconregs{i}(1:length(req.iconregs{i})-2) 2 1}};    
    end 
    wc=[]; 
end 
disp('Requested FL contrasts:'),  disp(req.cons(:,1))


% req.cons{1,2} =[  req.cons{1,2}; req.cons{2,2}]; disp('FAKE THING NOW'); 
            
            
if req.NonIdentityCons     
    error('Not done yet. set up in req.cons like the others.') 
    % req.cons{c,2} is a cell matrix where col 1= regname, col 1= weight on that regressor 


    for o1=1:1 % Load default specifications for requested model

        % Which contrast table to use?
        if sum(strcmp(log.onsetsmodel,{'m_c5_ChoiceRTFull_OULPEN';}))==1 % Competing + RT model
            contrasts.contrasttable='CompeteRT';

        else error('Error in Contrasts setup: Could not find requested onsets model. Which contrast table (i.e. excel sheet) to use?')
        end

        % Details for requested contrast table (excel file)
        switch contrasts.contrasttable
            case 'Compete'                    % Standard Choice + RL models ----------
                col.num_conds= 6+2*10 +4+4;
            case 'vMisc';
                col.num_conds=6+3+4+2+3;   disp('PredChoice #s of conds ok? If you added more, change settings');
            otherwise
                error('Error in Contrasts setup: Requested contrast table (i.e. excel sheet) not found')
        end
        col.contrastname=2;
        col.weightstart=3;
        row.conditionName=1;
        row.conditionType=2; % Condition or Pmod
        col.weightend=col.weightstart+col.num_conds-1;
        col.contrastnum=col.weightend+1;
        col.requested=col.contrastnum+1;

    end
    for o1=1:1 % Compile requested contrasts
    
    % Load default details regarding available contrasts
    [w.a1, w.a, w.req]=xlsread(req.contrasttable, contrasts.contrasttable);
    for i=col.contrastname+1:col.contrastnum-1 % Read details
        contrasts.condnames{i-col.contrastname}=w.req{row.conditionName,i};
        contrasts.condtypes{i-col.contrastname}=w.req{row.conditionType,i};
    end
    
    % Construct details for requested contrasts
    i=1;
    for r=row.conditionType+1:size(w.req,1);
        if isnan(w.req{r,col.requested})==0
            contrasts.contrastnames{i,1}=w.req{r,col.contrastname};
            contrasts.requested(i,1)=w.req{r,col.requested};
            
            % Load condition weights
            for c=col.contrastname+1:col.contrastnum-1 % Load all weights
                contrasts.weights(i,c-col.contrastname)=w.req{r,c};
            end
            %
            i=i+1;
        end
    end
      
     
end
end 
 

disp('Requested contrasts ##########################################')
disp(req.cons(:,1)) 
input('Continue to execute Contrasts?    ')

%% Execute Contrasts


for s=1:log.n_subjs 
    disp('Running contrasts #######################################')
    disp(['Subject ' num2str(s) '   (' log.subjects{s} ') -------------------'])
    ws.where=[where.data_brain filesep log.subjects{s} filesep '2 First level' log.FirstLevelThread filesep]; ws.wheremodel=[ws.where log.onsetsmodel ' Contrasted' filesep];
    if isdir(ws.wheremodel)==0; ws.estimfol=[ws.where log.onsetsmodel ' Estimated']; eval('java.io.File(ws.estimfol).renameTo(java.io.File(ws.wheremodel));');  else disp('Found contrasts folder. Assumed correct'); end
     f   = spm_select('List', ws.wheremodel, 'SPM.mat'); ws.spm  =load([ws.wheremodel f]); ws.regnamelist=ws.spm.SPM.xX.name';
    
    % Delete existing contrasts?
    if req.DeleteExistingFLContrasts && isfield(ws.spm.SPM, 'xCon')==1 && length(ws.spm.SPM.xCon)~=0; if s==1; input('Requested deletion of all existing FL contrasts. Continue?   '); end;  matlabbatch{1}.spm.stats.con.delete = 1;
    else; matlabbatch{1}.spm.stats.con.delete = 0;
    end
    
    % Specify requested contrasts
    disp('Specifying Contrasts ---------------------------');
    matlabbatch{1}.spm.stats.con.spmmat = cellstr([ws.wheremodel f]); c=1;
    for k= 1:size(req.cons,1) 
        matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';
        matlabbatch{1}.spm.stats.con.consess{c}.tcon.name = req.cons{k,1};
        %
        [ wc.regdetails ] = f_GetRegressorNums(ws.regnamelist, req.cons{k,2}); 
        wc.weights=zeros(1,length(ws.regnamelist)); 
        for j=1:size(req.cons{k,2},1);  % Find + weight requested regressors  
            wc.weights(wc.regdetails{j,2}) = req.cons{k,2}{j,3};  
        end 
        matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec=wc.weights; 
        wc=[]; c=c+1;
    end
    
    disp('Running contrasts -----------------')
    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    matlabbatch=[];wb=[];
    
end

%% END

disp('================c======================================'); w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
disp(' ');  disp('Analysis completed:')
disp(' '); disp(['No. of subjects: ' num2str(log.n_subjs)])
disp(' '); disp(['Model: ' log.onsetsmodel])
disp(' '); disp('Errors:'); disp(errorlog'); disp(' ')
disp('=======================================================')
 
diary off
try % Notify researcher
    f_sendemail('kurzlich', strcat(['PLPP Analysis batchscript is complete (s7_Firstlevel_2Contrasts  -  '  log.onsetsmodel  ')']), ' ',1);
end

 
%% Check if Contrasts are weighting the right regressors? 

do_check =0;
if do_check
    % Manually load SPM variable and execute
    for c= 1:length(SPM.xCon)
        disp(['Contrast:   ' SPM.xCon(c).name ' ------------------------'])
        disp([SPM.xX.name{find( SPM.xCon(c).c~=0 )}   num2str(SPM.xCon(c).c(SPM.xCon(c).c~=0 ))])
    end
end