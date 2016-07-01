function [ matlabbatch contrastfiles] = choice_cluster2x2(log,where,  subjectlist, firstlevelmodel, secondlevelmodel, Vars)
% Vars: {model name, {Vars{2}}, {FLcons and 2x2 assignment} 
% e.g.      Vars{1}: 'ConflictxChoice'
%             Vars{2}: {'cF' 'ct'; 'Rej' 'Expl'}  i.e. Fac1 Level1, Fac1 Level2, Fac2 Level1, Fac2 Level2,
%             Vars{3}: {'cF_Reject' [1 1];'cF_Explore'  [1 2];'ct_Bomb'  [2 1];'ct_Explore'  [2 2]};  

% You can manually assign things here too if you want to allocate.  

% Check inputs
if length(Vars)>3, error('Wrong no. inputs for SL! Request prob wrong'); end 


%% (1) Details for this model
%       Add to Vars{3}, col 3: name of contrast for that regressor (assume same for all subjects) 

f_mat=@(A,x)A(x);  
wc=load([where.data_brain filesep subjectlist{1} filesep '2 First level' log.FirstLevelThread filesep firstlevelmodel ' Contrasted' filesep 'SPM.mat']);
for c= 1:size(Vars{3},1) 
    Vars{3}{c,3}=  ['con_' f_mat(num2str(10000 + find(strcmp(cellstr(char(wc.SPM.xCon(:).name)), Vars{3}{c,1})) ), 2:5) '.img'];  
end
  
%% (2) Specify model for Factorial analysis
 
testfolder=  [where.resultsfolder filesep Vars{1} ]; 
if isdir(testfolder)==0; mkdir(testfolder); mkdir([testfolder filesep 'ROI']);end

% 2x2 ANOVA details in Vars{2}
matlabbatch{1}.spm.stats.factorial_design.dir = {testfolder};
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = [Vars{2}{1,1} ' v '  Vars{2}{1,2}];
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).dept = 1; % non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).variance = 0;% non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = [Vars{2}{2,1} ' v '  Vars{2}{2,2} ];
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).levels = 2;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).dept = 1;% non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).variance = 0;% non-default
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

% Specify Contrast files + Factorial cell 
contrastfiles=cell(size( Vars{3},1),1);
for i=1:size(Vars{3},1)
    matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).levels = Vars{3}{i,2}; % Levels in the factorial design 
    
    % Assign contrast files
    contrastfiles{i}=cell(length(subjectlist),1);
    for s=1:length(subjectlist)
        contrastfiles{i}{s}=[where.data_brain filesep subjectlist{s} filesep '2 First level' log.FirstLevelThread  filesep firstlevelmodel ' Contrasted' filesep  Vars{3}{i,3} ',1'];
    end
    matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).scans =contrastfiles{i};
end

% Include subjects as covariates
% matlabbatch{1}.spm.stats.factorial_design.cov(s).c = []; 
% matlabbatch{1}.spm.stats.factorial_design.cov(s).cname = [];
% matlabbatch{1}.spm.stats.factorial_design.cov(s).iCFI = 1;
% matlabbatch{1}.spm.stats.factorial_design.cov(s).iCC = 1;
% for s=1:length(subjectlist)
%     sub=zeros(1,length(subjectlist)); sub(s)=1;
%    matlabbatch{1}.spm.stats.factorial_design.cov(s).c = repmat(sub,[1 size(Vars{3},1)])';
%    matlabbatch{1}.spm.stats.factorial_design.cov(s).cname = ['sub_' subjectlist{s}];
%    matlabbatch{1}.spm.stats.factorial_design.cov(s).iCFI = 1;
%    matlabbatch{1}.spm.stats.factorial_design.cov(s).iCC = 1;
% end
% 
%% Run model

% Execute (Specify only)
disp('Specifying model --------------------------')
spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
matlabbatch=[];

% (2) Estimate model ##############################
disp('Estimating model ------------------------------')
matlabbatch{1}.spm.stats.fmri_est.spmmat = {[testfolder filesep 'SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
matlabbatch=[];

% (3) Add sensible contrasts
disp('Adding sensible contrasts ------------------------------')
matlabbatch{1}.spm.stats.con.spmmat = {[testfolder filesep 'SPM.mat']};
matlabbatch{1}.spm.stats.con.delete = 0; c=1;

% Identity matrix
matlabbatch{1}.spm.stats.con.consess{c}.fcon.name='Identity'; % Identity matrix (F Contrast)
matlabbatch{1}.spm.stats.con.consess{c}.fcon.convec ={[1 0 0 0]; [0 1 0 0];[0 0 1 0]; [0 0 0 1];};
matlabbatch{1}.spm.stats.con.consess{c}.fcon.sessrep = 'none'; c=c+1; 

% Sensible contrasts
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  [Vars{2}{1,1} ' > ' Vars{2}{1,2}]; % Factor 2
matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [1 1 -1 -1];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  [Vars{2}{1,1} ' < ' Vars{2}{1,2}];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [ -1 -1 1 1];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  [Vars{2}{2,1} ' > ' Vars{2}{2,2}]; % Factor 2 
matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [1 -1 1 -1];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name = [Vars{2}{2,1} ' < ' Vars{2}{2,2}];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [-1 1 -1 1];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;  
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  [Vars{2}{1,1} ' ' Vars{2}{2,1} ' > ' Vars{2}{2,2}]; % Fac 1 simple fx 
matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [1 -1 0 0];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  [Vars{2}{1,1} ' ' Vars{2}{2,1} ' < ' Vars{2}{2,2}];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [-1 1 0 0];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  [Vars{2}{1,2} ' ' Vars{2}{2,1} ' > ' Vars{2}{2,2}];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [0 0 1 -1 ];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  [Vars{2}{1,2} ' ' Vars{2}{2,1} ' < ' Vars{2}{2,2}];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [0 0 -1 1];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1; 
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  [Vars{2}{2,1} ' ' Vars{2}{1,1} ' > ' Vars{2}{1,2}]; % Fac 2 simple fx 
matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [1 0 -1 0];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  [Vars{2}{2,1} ' ' Vars{2}{1,1} ' < ' Vars{2}{1,2}];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [-1 0 1 0];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  [Vars{2}{2,2} ' ' Vars{2}{1,1} ' > ' Vars{2}{1,2}];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [0 1 0 -1];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;
matlabbatch{1}.spm.stats.con.consess{c}.tcon.name =  [Vars{2}{2,2} ' ' Vars{2}{1,1} ' < ' Vars{2}{1,2}];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.convec = [0 -1 0 1];
matlabbatch{1}.spm.stats.con.consess{c}.tcon.sessrep = 'none'; c=c+1;

%
spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
matlabbatch=[];

% (4) Save details to 2nd level folder
details.Vars=Vars;  details.confiles = contrastfiles;
save([testfolder filesep 'details_2ndlevel.mat'], 'details', 'log'); % Save details in 2nd level folder


end

