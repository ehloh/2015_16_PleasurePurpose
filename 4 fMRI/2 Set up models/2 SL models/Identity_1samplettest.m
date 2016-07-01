function [ matlabbatch contrastfiles] = Identity_1samplettest(log,where,  subjectlist, firstlevelmodel, secondlevelmodel, Vars,choices)
% Vars: {model/contrast names}   e.g. {'oe_event', 'trial'} ... 1sample ttest will be performed for all requested
 
%% (1) Details for this model

% Vars: Col 1=Name of model/contrast, Col 2: con file  
f_mat=@(A,x)A(x);  Vars = Vars(:); % Must be a column vector 
wc=load([where.data_brain filesep subjectlist{1} filesep '2 First level' log.FirstLevelThread filesep firstlevelmodel ' Contrasted' filesep 'SPM.mat']);
for v=1:size(Vars,1)
    if sum( strcmp(cellstr(char(wc.SPM.xCon(:).name)), Vars{v,1}) )~=1, error('No of found contrasts ~=1! Either too many or none were found.'), end
    Vars{v,2} = ['con_' f_mat(num2str(10000 + find( strcmp(cellstr(char(wc.SPM.xCon(:).name)), Vars{v,1})  )  ), 2:5) '.img'];  
end 
 
%% (3) Specify and estimate models for all Ttests

contrastfiles=cell(length(Vars),2);
for a=1:size(Vars,1)
        
        % (1) Specify model ##############################
        disp(['Specifying model no. '  num2str(a) '  :  ' Vars{a,1} '---------------------------'])
        testfolder=[where.resultsfolder filesep 'i ' Vars{a,1}];
        if isdir(testfolder)==0; mkdir(testfolder); mkdir([testfolder filesep 'ROI']); end

        matlabbatch{1}.spm.stats.factorial_design.dir = {testfolder};
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
        
        % Fetch subjects' scans/contrast files
        contrastfiles{a,1}=Vars{a,1};
        contrastfiles{a,2}=cell(length(subjectlist),1);
        for s=1:length(subjectlist)
            contrastfiles{a,2}{s}= [where.data_brain filesep subjectlist{s} filesep '2 First level' log.FirstLevelThread  filesep firstlevelmodel ' Contrasted' filesep Vars{a,2} ',1'];
        end
        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans=contrastfiles{a,2};

        % Execute (Specify only)
        spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
        matlabbatch=[];

        % (2) Estimate model ##############################
        disp(['Estimating model no. '  num2str(a) '  :  ' Vars{a,1} '---------------------------'])
        
        matlabbatch{1}.spm.stats.fmri_est.spmmat = {[testfolder filesep 'SPM.mat']};
        matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
        spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
        matlabbatch=[];
        
        
        % Save
        details.Vars=Vars(a,:);  details.confiles = contrastfiles; 
        save([testfolder filesep 'details_2ndlevel.mat'], 'details', 'log' ); % Save details in 2nd level folder

end
 
%% (4) Specify the contrast within each model 

disp('Specifying contrasts witin each model (only 1 available) ###############')

for  a=1:size(Vars,1)
       
    % Specify the contrast (only + and - available)
    matlabbatch{1}.spm.stats.con.spmmat = {[where.resultsfolder filesep 'i ' Vars{a,1} filesep 'SPM.mat']};
    matlabbatch{1}.spm.stats.con.delete = 0; 
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name =  'Pos';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = 1;
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name =  'Neg'; 
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec = -1;
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
    
    % Execute
    spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    matlabbatch=[];
    
end

end

