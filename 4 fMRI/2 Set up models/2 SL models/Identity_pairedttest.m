function [ matlabbatch contrastfiles] = Identity_pairedttest(log,where,  subjectlist, firstlevelmodel, secondlevelmodel, Vars,choices)
% Vars: {model/contrast name, {Cond1, Cond2}} 
%       e.g. {'oe_@@@rating', {'PL', 'PP'}: @@@ marker is replaced by PL & PP to search for the contrasts to compare w each other 

if length(Vars)~=2 , error('Wrong no. inputs for SL! Request prob wrong'); end 

%% (1) Details for this model 

% Instruc: Col 1=Name of comparison, Col 2&3: This minus That (names)
%               Col 4 & 5: This minus That (contrasts)
wc=load([where.data_brain filesep subjectlist{1} filesep '2 First level' log.FirstLevelThread filesep firstlevelmodel ' Contrasted' filesep 'SPM.mat']);
Instruc= cell(size(Vars,1), 5); f_mat=@(A,x)A(x);  
for v=1:size(size(Vars,1),1)
    Instruc{v,1} =  [ char(strrep(Vars{v,1}, '@@@', '')) ' ' Vars{v,2}{1} 'v' Vars{v,2}{2}];
    Instruc(v,2:3)= {strrep(Vars{v,1}, '@@@', 'PL' ) strrep(Vars{v,1}, '@@@', 'PP' )};
    
    % Find contrast numbers (assume constant across subjects)
    Instruc{v,4} =  ['con_' f_mat(num2str(10000 + find(strcmp(cellstr(char( wc.SPM.xCon.name)), Instruc{v,2})) ), 2:5) '.img']; 
    Instruc{v,5} =  ['con_' f_mat(num2str(10000 + find(strcmp(cellstr(char( wc.SPM.xCon.name)), Instruc{v,3})) ), 2:5) '.img'];  
end
 
  
%% (3) Specify and estimate models for all Ttests

contrastfiles=cell(size(Instruc,1),3);
for a=1:size(Instruc,1)
        
        % (1) Specify model ##############################
        disp(['Specifying model no. '  num2str(a) '  :  ' Instruc{a,1} '---------------------------'])
        testfolder=[where.resultsfolder filesep  Instruc{a,1}];
        if isdir(testfolder)==0; mkdir(testfolder); mkdir([testfolder filesep 'ROI']);end

        matlabbatch{1}.spm.stats.factorial_design.dir = {testfolder};
        matlabbatch{1}.spm.stats.factorial_design.des.pt.gmsca = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.pt.ancova = 0;
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
        
        % Fetch subjects' scans/contrast files
        contrastfiles{a,1}=Instruc{a,1};
        contrastfiles{a,2}=cell(length(subjectlist),1);
        contrastfiles{a,3}=cell(length(subjectlist),1);
        for s=1:length(subjectlist)
            contrastfiles{a,2}{s}= [where.data_brain filesep subjectlist{s} filesep '2 First level' log.FirstLevelThread filesep firstlevelmodel ' Contrasted' filesep Instruc{a,4} ',1'];
            contrastfiles{a,3}{s}= [where.data_brain filesep subjectlist{s} filesep '2 First level' log.FirstLevelThread filesep firstlevelmodel ' Contrasted' filesep Instruc{a,5} ',1'];
            matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(s).scans{1}=contrastfiles{a,2}{s};
            matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(s).scans{2}=contrastfiles{a,3}{s};
        end
        
        % Execute (Specify only)
        spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
        matlabbatch=[];

        % (2) Estimate model ##############################
        disp(['Estimating model no. '  num2str(a) '  :  ' Instruc{a,1} '---------------------------'])
        
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {[testfolder filesep 'SPM.mat']};
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
        matlabbatch=[];
        
        % Save details in 2nd level folder
        details.Vars=Vars;  details.confiles = contrastfiles; 
        save([testfolder filesep 'details_2ndlevel.mat'], 'details', 'log'); 

end


%% (4) Specify the contrast within each model 

disp('Specifying contrasts witin each model (only 1 available) ###############')

for  a=1:size(Instruc,1)
       
    % Weight subject covariates?
    wc=load([where.resultsfolder filesep  Instruc{a,1} filesep 'SPM.mat']);
    if size(wc.SPM.xX.X,2)==2+length(subjectlist), wc.subcovars=ones(1,length(subjectlist))*(1)/length(subjectlist);
    else wc.subcovars=[];
    end
            
    % Specify the contrast (only + and - available)
    matlabbatch{1}.spm.stats.con.spmmat = {[where.resultsfolder filesep Instruc{a,1} filesep 'SPM.mat']};
    matlabbatch{1}.spm.stats.con.delete = 0; 
    c=1;
    
    % START WITH WEIGHTS
    
         
    % (1) Overall
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name=[Vars{a,2}{1} ' and ' Vars{a,2}{2} 'pos'];
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [1 1 wc.subcovars*2];
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1; 
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name=[Vars{a,2}{1} ' and ' Vars{a,2}{2} 'neg'];
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [-1 -1 wc.subcovars*-2];
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;  
    
    % (2) Task-specific
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name=[Vars{a,2}{1} 'pos']; % First condition 
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [1 0 wc.subcovars*1];
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1; 
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name=[Vars{a,2}{1} 'neg'];
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [-1 0 wc.subcovars*-1];
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1; 
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name=[Vars{a,2}{2} 'pos']; % Second condition 
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [0 1 wc.subcovars*1];
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1; 
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name=[Vars{a,2}{2} 'neg'];
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [0 -1 wc.subcovars*-1];
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;  
    
    % (3) Difference beteen tasks 
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  [Vars{a,2}{1} ' > '  Vars{a,2}{2}]; % First > Second 
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [1 -1];
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1; 
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  [Vars{a,2}{1} ' < '  Vars{a,2}{2}]; % First < Second 
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [-1 1];
    matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none';c=c+1; 
    
    % Execute
    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    matlabbatch=[]; wc=[];
    
end

end

