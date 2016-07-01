% Preprocessing: realign & unwarp, coregister, segment, normalze, smooth
clear all;close all hidden; clc

% Requested analysis
request.realignunwarp=0;
request.coregister=0;
request.normalize=1;
request.smooth=1;  
request.smoothingsize=8; % mm
log.specificsubjects={     
%     'p01_YH';  'p02_MI'; 
% 'p03_AY';
% 'p04_AW'; 'p05_CA';'p06_BT';'p07_HC';'p08_KC'; 'p09_KJ'; 'p10_YC'; 
%     'p11_BI'; 
'p12_AL'; 
% 'p13_MS';'p14_SK';'p15_MM';'p16_SH';'p17_BB';'p18_WB'; 'p19_HB';'p20_LZ'
}; % BLANK to request all subjects

% Remember to turn on try catches 
for o1=1:1 % General settings and specifications
    
    % Load subjects
    w=pwd;  if strcmp(w(1), '/')==1;  where.where='/Users/EleanorL/Dropbox/SCRIPPS/3b PLPP fmri';  where.data_brain='/Users/EleanorL/Desktop/3 PLPR';  
    
    where.spm='/Users/EleanorL/Documents/MATLAB/spm8';
    else where.where='C:\Users\e.loh\Dropbox\SCRIPPS\3b PLPP fmri';  where.data_brain= 'D:\1 PLPP\1 MRI data'; 
        where.spm='C:\toolbox\64\spm8';
    end   
    [n t log.datalog]=xlsread([where.data_brain fs  'datalog_plpr.xlsx']); addpath(where.where)
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    col.funcsess=[7 8 9];  col.funcfieldmaps=[10 11 12];  col.struct=13;  col.func_nscans=[14 15 16];  % log.datalog contents 
    scan=load([where.where filesep '1 Preprocessing' fs 'i_scanningdetails.mat']);
    errorlog=cell(1,1); e=1;  % Log
    
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

%% Step 1: Realign & unwarp (prefix u)

if request.realignunwarp==1
    disp(' ########## (1A) REALIGN & UNWARP: Create VDM files ##########')
    default.fieldmap=[where.spm filesep 'toolbox' filesep 'FieldMap' filesep 'pm_defaults_Sonata_eFoV.m'];  % Checked w Caroline 
    do_vdm=1;
    if do_vdm
        for s=1:length(log.subjects)  % Create vdm file
            disp(['Subject ' num2str(s)   '   (' log.subjects{s} ')  --------------- '])
            if strcmp(log.subjects{s}, 'p12_AL')==1, disp('p12_AL skipped: no fieldmap'); else 
                try
                    ws.where=[where.data_brain filesep log.subjects{s} filesep '1 Preprocessed' filesep];
                    % Check: do we have VDMs already? If so, omit. 
                    ws.scid= spm_select('List', [ws.where 'Func_s1'], '^bf*' );  ws.scid= ws.scid(1,1:strfind(ws.scid(1,:), '-'));
                    scanid{s,1}=ws.scid; ws.fid=  ws.scid(strfind(ws.scid,'fB'):end);  f= [spm_select('List',[ws.where 'Fieldmap'], ['^vdm5_scs' ws.fid(2:end) '.*session1.img']);spm_select('List',[ws.where 'Fieldmap'], ['^vdm5_scs' ws.fid(2:end) '.*session2.img']); spm_select('List',[ws.where 'Fieldmap'], ['^vdm5_scs' ws.fid(2:end) '.*session3.img']);spm_select('List',[ws.where 'Fieldmap'], ['^vdm5_scs' ws.fid(2:end) '.*session4.img']); spm_select('List',[ws.where 'Fieldmap'], ['^vdm5_scs' ws.fid(2:end) '.*session5.img']); spm_select('List',[ws.where 'Fieldmap'], ['^vdm5_scs' ws.fid(2:end) '.*session6.img']);];
                    if size(f,1)==5; errorlog{e,1}=['Log non-failure: (1a) VDMs already exist, so not created --- ' log.datalog{s+1,1} ]; disp(errorlog{e,1}); e=e+1;
                    else
                        
                        % Create VDM files
                        f=spm_select('List',[ws.where 'Fieldmap' fs],'.*img$');  % 1=mag, 2=phase
                        matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.phase=cellstr([ws.where 'Fieldmap' fs strtrim( f(1,:) ) ',1']);
                        matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.magnitude=cellstr([ws.where 'Fieldmap' fs strtrim( f(2,:) ) ',1']);
                        
                        % first epi for each run
                        ws.epis=  [repmat([ws.where 'Func_s1' filesep],  size(spm_select('List',[ws.where 'Func_s1' filesep], ['^bf.*00' num2str(scan.nDummyVols+1) '-01.img$']),1), 1) spm_select('List',[ws.where 'Func_s1' filesep], ['^bf.*00' num2str(scan.nDummyVols+1) '-01.img$'])]  ;
                        matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.session(1).epi=cellstr(ws.epis(1,:)); 
                        ws.epis=  [repmat([ws.where 'Func_s2' filesep],  size(spm_select('List',[ws.where 'Func_s2' filesep], ['^bf.*00' num2str(scan.nDummyVols+1) '-01.img$']),1), 1) spm_select('List',[ws.where 'Func_s2' filesep], ['^bf.*00' num2str(scan.nDummyVols+1) '-01.img$'])]  ;
                        matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.session(2).epi=cellstr(ws.epis(1,:)); 
                        ws.epis=  [repmat([ws.where 'Func_s3' filesep],  size(spm_select('List',[ws.where 'Func_s3' filesep], ['^bf.*00' num2str(scan.nDummyVols+1) '-01.img$']),1), 1) spm_select('List',[ws.where 'Func_s3' filesep], ['^bf.*00' num2str(scan.nDummyVols+1) '-01.img$'])]  ;
                        matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.session(3).epi=cellstr(ws.epis(1,:)); 
                        matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.defaults.defaultsfile = {default.fieldmap}; 
                        matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.matchvdm = 1;
                        matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.sessname = 'session';
                        matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.writeunwarped = 0;
                        matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.anat = '';
                        matlabbatch{1}.spm.tools.fieldmap.presubphasemag.subj.matchanat = 0;
%                         save('batch', 'matlabbatch') 
                        spm_jobman('initcfg'); spm_jobman('run',matlabbatch); matlabbatch=[];  % Run ! 
                    end
                catch, errorlog{e,1}=['Failed: (1a) Realign & Unwarp: Create VDM file  --- ' log.datalog{s+1,1} ]; disp(errorlog{e,1}); e=e+1;
                end
            end
            ws=[];
        end
        
        disp(' ########## (1B) REALIGN & UNWARP: Execute realigning & unwarping ##########')
        for o1 =1:1 % Realign & unwarp settings
            settings.realignunwarp.eoptions.quality = 0.9;
            settings.realignunwarp.eoptions.sep = 4;
            settings.realignunwarp.eoptions.fwhm = 5;
            settings.realignunwarp.eoptions.rtm = 0;
            settings.realignunwarp.eoptions.einterp = 2;
            settings.realignunwarp.eoptions.ewrap = [0 0 0];
            settings.realignunwarp.eoptions.weight = '';
            settings.realignunwarp.uweoptions.basfcn = [12 12];
            settings.realignunwarp.uweoptions.regorder = 1;
            settings.realignunwarp.uweoptions.lambda = 100000;
            settings.realignunwarp.uweoptions.jm = 0;
            settings.realignunwarp.uweoptions.fot = [4 5];
            settings.realignunwarp.uweoptions.sot = [];
            settings.realignunwarp.uweoptions.uwfwhm = 4;
            settings.realignunwarp.uweoptions.rem = 1;
            settings.realignunwarp.uweoptions.noi = 5;
            settings.realignunwarp.uweoptions.expround = 'Average';
            settings.realignunwarp.uwroptions.uwwhich = [2 1];
            settings.realignunwarp.uwroptions.rinterp = 4;
            settings.realignunwarp.uwroptions.wrap = [0 0 0];
            settings.realignunwarp.uwroptions.mask = 1;
            settings.realignunwarp.uwroptions.uwroptions.prefix = 'u';
        end
        for s=1:length(log.subjects)
            disp(['Subject ' num2str(s) '  (' log.subjects{s} ') --------------- '])
            if strcmp(log.subjects{s}, 'p12_AL')==1, disp('p12_AL skipped: no fieldmap'); else
                try

                    ws.where=[where.data_brain filesep log.datalog{s+1,1} filesep '1 Preprocessed' filesep];
                    ws.scid= spm_select('List', [ws.where 'Func_s1'], '^bf*' );  ws.scid= ws.scid(1,1:strfind(ws.scid(1,:), '-'));
                    scanid{s,1}=ws.scid; ws.fid=  ws.scid(strfind(ws.scid,'fB'):end);  
                    
                    % For each run, select EPIs & vdm file (from fieldmap folder)
                    matlabbatch{1}.spm.spatial.realignunwarp=settings.realignunwarp;
                    ss=1;k=1;   ws.epis= num2str(10000+log.datalog{s+1, col.funcsess(ss)});   ws.epis=  spm_select('List', [ws.where 'Func_s' num2str(ss)], ['^' ws.scid ws.epis(2:end)  '.*.img' ]);
                    matlabbatch{1}.spm.spatial.realignunwarp.data(k).scans = cellstr([repmat([ws.where 'Func_s' num2str(ss) fs],  size(ws.epis,1),1) ws.epis repmat(',1',  size(ws.epis,1),1)]);
                    matlabbatch{1}.spm.spatial.realignunwarp.data(k).pmscan = { [ws.where 'Fieldmap'  fs  spm_select('List', [ws.where 'Fieldmap'], ['^vdm5_scs' ws.fid(2:end) '.*session' num2str(k) '.img'])  ',1']};  
                    ss=2;k=k+1;  ws.epis= num2str(10000+log.datalog{s+1, col.funcsess(ss)});   ws.epis=  spm_select('List', [ws.where 'Func_s' num2str(ss)], ['^' ws.scid ws.epis(2:end)  '.*.img' ]);
                    matlabbatch{1}.spm.spatial.realignunwarp.data(k).scans = cellstr([repmat([ws.where 'Func_s' num2str(ss) fs],  size(ws.epis,1),1) ws.epis repmat(',1',  size(ws.epis,1),1)]);
                    matlabbatch{1}.spm.spatial.realignunwarp.data(k).pmscan = { [ws.where 'Fieldmap'  fs  spm_select('List', [ws.where 'Fieldmap'], ['^vdm5_scs' ws.fid(2:end) '.*session' num2str(k) '.img'])  ',1']};  
                    ss=3;k=k+1; ws.epis= num2str(10000+log.datalog{s+1, col.funcsess(ss)});   ws.epis=  spm_select('List', [ws.where 'Func_s' num2str(ss)], ['^' ws.scid ws.epis(2:end)  '.*.img' ]);
                    matlabbatch{1}.spm.spatial.realignunwarp.data(k).scans = cellstr([repmat([ws.where 'Func_s' num2str(ss) fs],  size(ws.epis,1),1) ws.epis repmat(',1',  size(ws.epis,1),1)]);
                    matlabbatch{1}.spm.spatial.realignunwarp.data(k).pmscan = { [ws.where 'Fieldmap'  fs  spm_select('List', [ws.where 'Fieldmap'], ['^vdm5_scs' ws.fid(2:end) '.*session' num2str(k) '.img'])  ',1']};   
                    spm_jobman('run' , matlabbatch);                     matlabbatch=[];                    
                catch; errorlog{e,1}=['Failed: (1b) Realign & Unwarp  --- ' log.datalog{s+1,1}]; disp(errorlog{e,1}); e=e+1;
                end
            end
            ws=[];
        end
    end
    
    
    % Realign only for p12 (prefix = r)
    if strcmp(log.subjects,'p12_AL')
        disp('Manual realignment for p12 ------------------------------------------------') 
        for o=1:1  % Manual collection of files for session 1 and 2
            
            matlabbatch{1}.spm.spatial.realign.estwrite.data = {
                {
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00007-000007-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00008-000008-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00009-000009-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00010-000010-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00011-000011-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00012-000012-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00013-000013-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00014-000014-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00015-000015-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00016-000016-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00017-000017-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00018-000018-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00019-000019-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00020-000020-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00021-000021-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00022-000022-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00023-000023-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00024-000024-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00025-000025-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00026-000026-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00027-000027-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00028-000028-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00029-000029-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00030-000030-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00031-000031-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00032-000032-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00033-000033-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00034-000034-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00035-000035-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00036-000036-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00037-000037-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00038-000038-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00039-000039-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00040-000040-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00041-000041-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00042-000042-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00043-000043-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00044-000044-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00045-000045-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00046-000046-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00047-000047-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00048-000048-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00049-000049-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00050-000050-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00051-000051-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00052-000052-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00053-000053-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00054-000054-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00055-000055-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00056-000056-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00057-000057-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00058-000058-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00059-000059-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00060-000060-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00061-000061-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00062-000062-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00063-000063-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00064-000064-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00065-000065-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00066-000066-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00067-000067-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00068-000068-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00069-000069-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00070-000070-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00071-000071-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00072-000072-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00073-000073-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00074-000074-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00075-000075-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00076-000076-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00077-000077-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00078-000078-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00079-000079-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00080-000080-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00081-000081-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00082-000082-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00083-000083-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00084-000084-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00085-000085-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00086-000086-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00087-000087-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00088-000088-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00089-000089-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00090-000090-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00091-000091-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00092-000092-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00093-000093-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00094-000094-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00095-000095-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00096-000096-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00097-000097-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00098-000098-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00099-000099-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00100-000100-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00101-000101-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00102-000102-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00103-000103-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00104-000104-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00105-000105-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00106-000106-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00107-000107-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00108-000108-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00109-000109-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00110-000110-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00111-000111-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00112-000112-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00113-000113-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00114-000114-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00115-000115-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00116-000116-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00117-000117-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00118-000118-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00119-000119-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00120-000120-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00121-000121-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00122-000122-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00123-000123-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00124-000124-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00125-000125-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00126-000126-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00127-000127-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00128-000128-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00129-000129-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00130-000130-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00131-000131-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00132-000132-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00133-000133-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00134-000134-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00135-000135-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00136-000136-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00137-000137-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00138-000138-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00139-000139-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00140-000140-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00141-000141-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00142-000142-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00143-000143-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00144-000144-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00145-000145-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00146-000146-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00147-000147-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00148-000148-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00149-000149-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00150-000150-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00151-000151-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00152-000152-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00153-000153-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00154-000154-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00155-000155-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00156-000156-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00157-000157-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00158-000158-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00159-000159-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00160-000160-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00161-000161-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00162-000162-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00163-000163-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00164-000164-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00165-000165-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00166-000166-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00167-000167-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00168-000168-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00169-000169-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00170-000170-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00171-000171-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00172-000172-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00173-000173-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00174-000174-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00175-000175-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00176-000176-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00177-000177-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00178-000178-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00179-000179-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00180-000180-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00181-000181-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00182-000182-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00183-000183-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00184-000184-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00185-000185-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00186-000186-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00187-000187-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00188-000188-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00189-000189-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00190-000190-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00191-000191-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00192-000192-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00193-000193-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00194-000194-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00195-000195-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00196-000196-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00197-000197-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00198-000198-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00199-000199-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00200-000200-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00201-000201-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00202-000202-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00203-000203-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00204-000204-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00205-000205-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00206-000206-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00207-000207-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00208-000208-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00209-000209-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00210-000210-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00211-000211-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00212-000212-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00213-000213-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00214-000214-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00215-000215-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00216-000216-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00217-000217-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00218-000218-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00219-000219-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00220-000220-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00221-000221-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00222-000222-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00223-000223-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00224-000224-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00225-000225-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00226-000226-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00227-000227-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00228-000228-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00229-000229-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00230-000230-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00231-000231-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00232-000232-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00233-000233-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00234-000234-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00235-000235-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00236-000236-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00237-000237-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00238-000238-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00239-000239-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00240-000240-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00241-000241-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00242-000242-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00243-000243-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00244-000244-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00245-000245-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00246-000246-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00247-000247-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00248-000248-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00249-000249-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00250-000250-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00251-000251-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00252-000252-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00253-000253-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00254-000254-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00255-000255-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00256-000256-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00257-000257-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00258-000258-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00259-000259-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00260-000260-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00261-000261-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00262-000262-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00263-000263-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00264-000264-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00265-000265-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00266-000266-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00267-000267-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00268-000268-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00269-000269-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00270-000270-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00271-000271-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00272-000272-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00273-000273-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00274-000274-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00275-000275-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00276-000276-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00277-000277-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00278-000278-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00279-000279-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00280-000280-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00281-000281-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00282-000282-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00283-000283-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00284-000284-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00285-000285-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00286-000286-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00287-000287-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00288-000288-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00289-000289-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00290-000290-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00291-000291-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00292-000292-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00293-000293-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00294-000294-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00295-000295-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00296-000296-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00297-000297-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00298-000298-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00299-000299-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00300-000300-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00301-000301-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00302-000302-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00303-000303-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00304-000304-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00305-000305-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00306-000306-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00307-000307-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00308-000308-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00309-000309-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00310-000310-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00311-000311-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00312-000312-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00313-000313-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00314-000314-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00315-000315-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00316-000316-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00317-000317-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00318-000318-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00319-000319-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00320-000320-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00321-000321-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00322-000322-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00323-000323-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00324-000324-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00325-000325-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00326-000326-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00327-000327-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00328-000328-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00329-000329-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00330-000330-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00331-000331-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00332-000332-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00333-000333-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00334-000334-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00335-000335-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00336-000336-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00337-000337-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00338-000338-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00339-000339-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00340-000340-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00341-000341-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00342-000342-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00343-000343-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00344-000344-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00345-000345-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00346-000346-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00347-000347-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00348-000348-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00349-000349-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00350-000350-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00351-000351-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00352-000352-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00353-000353-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00354-000354-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00355-000355-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00356-000356-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00357-000357-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s1\bfB005200-0003-00358-000358-01.img,1'
                }
                {
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00007-000007-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00008-000008-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00009-000009-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00010-000010-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00011-000011-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00012-000012-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00013-000013-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00014-000014-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00015-000015-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00016-000016-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00017-000017-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00018-000018-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00019-000019-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00020-000020-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00021-000021-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00022-000022-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00023-000023-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00024-000024-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00025-000025-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00026-000026-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00027-000027-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00028-000028-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00029-000029-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00030-000030-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00031-000031-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00032-000032-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00033-000033-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00034-000034-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00035-000035-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00036-000036-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00037-000037-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00038-000038-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00039-000039-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00040-000040-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00041-000041-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00042-000042-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00043-000043-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00044-000044-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00045-000045-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00046-000046-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00047-000047-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00048-000048-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00049-000049-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00050-000050-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00051-000051-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00052-000052-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00053-000053-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00054-000054-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00055-000055-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00056-000056-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00057-000057-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00058-000058-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00059-000059-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00060-000060-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00061-000061-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00062-000062-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00063-000063-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00064-000064-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00065-000065-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00066-000066-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00067-000067-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00068-000068-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00069-000069-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00070-000070-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00071-000071-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00072-000072-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00073-000073-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00074-000074-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00075-000075-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00076-000076-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00077-000077-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00078-000078-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00079-000079-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00080-000080-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00081-000081-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00082-000082-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00083-000083-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00084-000084-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00085-000085-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00086-000086-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00087-000087-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00088-000088-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00089-000089-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00090-000090-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00091-000091-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00092-000092-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00093-000093-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00094-000094-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00095-000095-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00096-000096-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00097-000097-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00098-000098-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00099-000099-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00100-000100-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00101-000101-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00102-000102-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00103-000103-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00104-000104-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00105-000105-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00106-000106-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00107-000107-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00108-000108-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00109-000109-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00110-000110-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00111-000111-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00112-000112-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00113-000113-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00114-000114-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00115-000115-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00116-000116-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00117-000117-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00118-000118-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00119-000119-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00120-000120-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00121-000121-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00122-000122-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00123-000123-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00124-000124-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00125-000125-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00126-000126-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00127-000127-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00128-000128-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00129-000129-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00130-000130-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00131-000131-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00132-000132-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00133-000133-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00134-000134-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00135-000135-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00136-000136-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00137-000137-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00138-000138-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00139-000139-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00140-000140-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00141-000141-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00142-000142-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00143-000143-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00144-000144-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00145-000145-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00146-000146-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00147-000147-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00148-000148-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00149-000149-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00150-000150-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00151-000151-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00152-000152-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00153-000153-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00154-000154-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00155-000155-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00156-000156-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00157-000157-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00158-000158-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00159-000159-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00160-000160-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00161-000161-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00162-000162-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00163-000163-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00164-000164-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00165-000165-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00166-000166-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00167-000167-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00168-000168-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00169-000169-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00170-000170-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00171-000171-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00172-000172-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00173-000173-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00174-000174-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00175-000175-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00176-000176-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00177-000177-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00178-000178-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00179-000179-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00180-000180-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00181-000181-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00182-000182-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00183-000183-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00184-000184-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00185-000185-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00186-000186-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00187-000187-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00188-000188-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00189-000189-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00190-000190-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00191-000191-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00192-000192-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00193-000193-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00194-000194-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00195-000195-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00196-000196-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00197-000197-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00198-000198-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00199-000199-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00200-000200-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00201-000201-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00202-000202-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00203-000203-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00204-000204-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00205-000205-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00206-000206-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00207-000207-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00208-000208-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00209-000209-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00210-000210-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00211-000211-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00212-000212-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00213-000213-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00214-000214-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00215-000215-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00216-000216-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00217-000217-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00218-000218-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00219-000219-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00220-000220-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00221-000221-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00222-000222-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00223-000223-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00224-000224-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00225-000225-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00226-000226-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00227-000227-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00228-000228-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00229-000229-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00230-000230-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00231-000231-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00232-000232-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00233-000233-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00234-000234-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00235-000235-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00236-000236-01.img,1'
                'D:\1 PLPP\1 MRI data\p12_AL\1 Preprocessed\Func_s2\bfB005200-0005-00237-000237-01.img,1'
                }
                }';
            
        end
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 4;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
        matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 2;
        spm_jobman('run' , matlabbatch); matlabbatch=[];                    
    end
    
end

%% Step 2: Coregister (no prefix)

if request.coregister==1
    disp(' ############### (2) COREGISTRATION ############ ##########')
    for o2=1:1 % Settings for Coregistration
        settings.coreg.other = {''};
        settings.coreg.eoptions.cost_fun = 'nmi';
        settings.coreg.eoptions.sep = [4 2];
        settings.coreg.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
        settings.coreg.eoptions.fwhm = [7 7];
    end
    for s=1:length(log.subjects)
        try
            disp(['Subject ' num2str(s) '  (' log.subjects{s} ') --------------- '])
            matlabbatch{1}.spm.spatial.coreg.estimate=settings.coreg;  % Specifications 
            ws.where=[where.data_brain filesep log.datalog{s+1,1} filesep '1 Preprocessed' filesep];
            if strcmp(log.subjects{s}, 'p12_AL')==1,  ws.pref='^rbf*' ;  else ws.pref='^ubf*' ; end
            ws.scid= spm_select('List', [ws.where 'Func_s1'], ws.pref);   ws.scid= ws.scid(1,1:strfind(ws.scid(1,:), '-'));
                
            % Choose files
            matlabbatch{1}.spm.spatial.coreg.estimate.ref={[ws.where 'Structural' filesep spm_select('List', [ws.where 'Structural'], '^s*.*.img') ',1']};  % Reference: Structural
            matlabbatch{1}.spm.spatial.coreg.estimate.source={[ws.where 'Func_s1' fs spm_select('List', [ws.where 'Func_s1'], ['^' ws.scid f_subindex(num2str(10000 + log.datalog{s+1, col.funcsess(1)}), 2:5) '.*.' f_subindex(num2str(1000000+ scan.nDummyVols+1), 2:7) '-01.img' ]) ',1']};   % Source image: meanubfM file (or, 1st volume from 1st run) 
            matlabbatch{1}.spm.spatial.coreg.estimate.other=[  % Collect functionals from all runs
                [repmat([ws.where 'Func_s1' fs ], size(spm_select('List', [ws.where 'Func_s1'],  ['^' ws.scid '*.*img']),1),1)  spm_select('List', [ws.where 'Func_s1'],  ['^' ws.scid '*.*img'])]
                [repmat([ws.where 'Func_s2' fs ], size(spm_select('List', [ws.where 'Func_s2'],  ['^' ws.scid '*.*img']),1),1)  spm_select('List', [ws.where 'Func_s2'],  ['^' ws.scid '*.*img'])]
                [repmat([ws.where 'Func_s3' fs ], size(spm_select('List', [ws.where 'Func_s3'],  ['^' ws.scid '*.*img']),1),1)  spm_select('List', [ws.where 'Func_s3'],  ['^' ws.scid '*.*img'])]
                ]; 
            matlabbatch{1}.spm.spatial.coreg.estimate.other= cellstr([matlabbatch{1}.spm.spatial.coreg.estimate.other repmat(',1', size(matlabbatch{1}.spm.spatial.coreg.estimate.other,1),1)]); 
            spm_jobman('run' , matlabbatch); matlabbatch=[];
        catch; errorlog{e,1}=['Failed: (2) Coregister & Reslice --- ' log.datalog{s+1,1}]; e=e+1;
        end
        ws=[];
    end
end

%% Step 3: Normalization (prefix 'w')
% Normalization here does NOT use the MPMs (which involves segmentation)
%   See other scripts for MPM-normalization.

if request.normalize==1  
    disp(' ############### (3) Normalization: Execute normalization ############ ##########')
    for o1=1:1 % Settings for Normalization 
        settings.normalization.roptions.prefix = 'w';
        settings.normalization.eoptions.weight = '';
        settings.normalization.eoptions.smosrc = 8;  % Source image smoothing; should match the template (SPM templates are smoothed to 8mm)
        settings.normalization.eoptions.smoref = 0; % Additional smoothing for template (standard space)
        settings.normalization.eoptions.regtype = 'mni';
        settings.normalization.eoptions.cutoff = 25;
        settings.normalization.eoptions.nits = 16;
        settings.normalization.eoptions.reg = 1;
        settings.normalization.roptions.preserve = 0;
        settings.normalization.roptions.bb = [-78 -112 -50;  78 76 85];
        settings.normalization.roptions.interp = 7; 
        settings.normalization.roptions.vox = [3 3 2]; % Effective voxel sizes (check w physics) 
        settings.normalization.eoptions.template = {[where.spm fs 'templates' fs 'EPI.nii,1']};
        settings.normalization.roptions.wrap = [0 1 0];        
    end
    for s=1:length(log.subjects)
%         try
            disp(['Subject ' num2str(s) '  (' log.subjects{s} ') --------------- '])
            ws.where=[where.data_brain fs log.subjects{s} fs  '1 Preprocessed' filesep];
            if strcmp(log.subjects{s}, 'p12_AL')==1,  ws.pref='^rbf*' ;  else ws.pref='^ubf*' ; end
            ws.scid= spm_select('List', [ws.where 'Func_s1'], ws.pref);   ws.scid= ws.scid(1,1:strfind(ws.scid(1,:), '-'));
            matlabbatch{1}.spm.spatial.normalise.estwrite=settings.normalization;
            matlabbatch{1}.spm.spatial.normalise.estwrite.subj.wtsrc = '';
            
%            
%             [repmat([ws.where 'Func_s1' fs ], size(spm_select('List', [ws.where 'Func_s1'],  ['^' ws.scid '*.*img']),1),1)  spm_select('List', [ws.where 'Func_s1'],  ['^' ws.scid '*.*img'])]
%                 [repmat([ws.where 'Func_s2' fs ], size(spm_select('List', [ws.where 'Func_s2'],  ['^' ws.scid '*.*img']),1),1)  spm_select('List', [ws.where 'Func_s2'],  ['^' ws.scid '*.*img'])]
%                 [repmat([ws.where 'Func_s3' fs ], size(spm_select('List', [ws.where 'Func_s3'],  ['^' ws.scid '*.*img']),1),1)  spm_select('List', [ws.where 'Func_s3'],  ['^' ws.scid '*.*img'])]
               
            
            
            % Load scans
            matlabbatch{1}.spm.spatial.normalise.estwrite.subj.source =   {[ws.where 'Func_s1' fs spm_select('List', [ws.where 'Func_s1'], ['^' ws.scid f_subindex(num2str(10000 + log.datalog{s+1, col.funcsess(1)}), 2:5) '.*.' f_subindex(num2str(1000000+ scan.nDummyVols+1), 2:7) '-01.img' ]) ',1']};   % Source image: meanubfM file (or, 1st volume from 1st run)
            matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = [  % Collect functionals from all runs
                [repmat([ws.where 'Func_s1' fs ], size(spm_select('List', [ws.where 'Func_s1'],  ['^' ws.scid '*.*img']),1),1)  spm_select('List', [ws.where 'Func_s1'],  ['^' ws.scid '*.*img'])]
                [repmat([ws.where 'Func_s2' fs ], size(spm_select('List', [ws.where 'Func_s2'],  ['^' ws.scid '*.*img']),1),1)  spm_select('List', [ws.where 'Func_s2'],  ['^' ws.scid '*.*img'])]
                [repmat([ws.where 'Func_s3' fs ], size(spm_select('List', [ws.where 'Func_s3'],  ['^' ws.scid '*.*img']),1),1)  spm_select('List', [ws.where 'Func_s3'],  ['^' ws.scid '*.*img'])]
                ];
            matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = cellstr([matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample  repmat(',1', size(matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample ,1),1)]);
            spm_jobman('run' , matlabbatch); matlabbatch=[];
%         catch; errorlog{e,1}=['Failed: (3) Normalization  --- ' log.datalog{s+1,1}]; e=e+1;
%         end
        ws=[];
    end
end 

%% Step 5: Smoothing
% Always saved in different folders

if request.smooth==1  
    disp([' ############### (4) SMOOTHING (size: ' num2str(request.smoothingsize) 'x' num2str(request.smoothingsize) 'x' num2str(request.smoothingsize) ')   ############ ##########'])
    for o1=1:1 % Settings for Smoothing
        settings.smooth.fwhm = [request.smoothingsize request.smoothingsize request.smoothingsize]; % Smoothing size
        settings.smooth.dtype = 0;
        settings.smooth.im = 0;
        settings.smooth.prefix = ['s' num2str(request.smoothingsize)];
    end
    for s=1:length(log.subjects)
        try 
            disp(['Subject ' num2str(s) '  (' log.subjects{s} ') --------------- '])
            matlabbatch{1}.spm.spatial.smooth=settings.smooth;
            ws.where=[where.data_brain fs log.subjects{s} fs '1 Preprocessed' fs];
            if strcmp(log.subjects{s}, 'p12_AL')==1,  ws.pref='^wrbf*' ;  else ws.pref='^wubf*' ; end
            ws.scid= spm_select('List', [ws.where 'Func_s1'], ws.pref);   ws.scid= ws.scid(1,1:strfind(ws.scid(1,:), '-')); 
            
            matlabbatch{1}.spm.spatial.smooth.data = [  % Collect functionals from all runs
                [repmat([ws.where 'Func_s1' fs ], size(spm_select('List', [ws.where 'Func_s1'],  ['^' ws.scid '*.*img']),1),1)  spm_select('List', [ws.where 'Func_s1'],  ['^' ws.scid '*.*img'])]
                [repmat([ws.where 'Func_s2' fs ], size(spm_select('List', [ws.where 'Func_s2'],  ['^' ws.scid '*.*img']),1),1)  spm_select('List', [ws.where 'Func_s2'],  ['^' ws.scid '*.*img'])]
                [repmat([ws.where 'Func_s3' fs ], size(spm_select('List', [ws.where 'Func_s3'],  ['^' ws.scid '*.*img']),1),1)  spm_select('List', [ws.where 'Func_s3'],  ['^' ws.scid '*.*img'])]
                ];
            matlabbatch{1}.spm.spatial.smooth.data= cellstr([ matlabbatch{1}.spm.spatial.smooth.data repmat(',1', size(matlabbatch{1}.spm.spatial.smooth.data ,1),1)]); 
            spm_jobman('run' , matlabbatch); matlabbatch=[];
        catch; errorlog{e,1}=['Failed: (4) Smoothing (size: ' num2str(request.smoothingsize) ' --- ' log.datalog{s+1,1}]; e=e+1;
        end
        ws=[];
    end
end
 

%% END

disp('=======================================================')
w.c=clock; disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
disp(' '); disp('Analysis completed:')
disp(request); disp(['No. of subjects: ' num2str(log.n_subjs)])
disp(' '); disp(errorlog); disp(' ')
disp('=======================================================')
try  f_sendemail('kurzlich', strcat('Analysis batchscript is complete (s3_Preprocessing)'), ' ',1); end


% disp('Log scanids in datalog file!'), disp([log.subjects scanid])