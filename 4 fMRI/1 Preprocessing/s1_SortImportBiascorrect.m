% SortImportDeletedummiesBiascorrect
clear all;close all hidden; clc

% Requested analysis
request.sort=1;
request.deletedummyvolumes=1;
request.biascorrect=1;
%   
log.specificsubjects={ 
    'p12_AL'
    
% 'p01_YH';  'p02_MI';
% 'p03_AY'; 'p04_AW';'p05_CA';'p06_BT';'p07_HC';'p08_KC';'p09_KJ';   'p10_YC'; 'p11_BI';  'p12_AL';'p13_MS';'p14_SK';'p15_MM';'p16_SH';'p17_BB';'p18_WB';'p19_HB';'p20_LZ'

}; % BLANK to request all subjects

for o1=1:1 % General settings and specifications
    
     % Load subjects
    w=pwd;  if strcmp(w(1), '/')==1;  where.where='/Users/EleanorL/Dropbox/SCRIPPS/3b PLPP fmri';  where.data_brain='/Users/EleanorL/Desktop/3 PLPR';  
    else where.where='C:\Users\e.loh\Dropbox\SCRIPPS\3b PLPP fmri';  where.data_brain= 'D:\1 PLPP\1 MRI data'; 
    end  
    addpath(where.where)
    [n t log.datalog]=xlsread([where.data_brain filesep 'datalog_plpr.xlsx']);
%     [n t log.datalog]=xlsread('/Users/EleanorL/Dropbox/SCRIPPS/7 Pleasure purpose/3 Behaviour/datalog_plpr.xlsx');
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    col.funcsess=[7 8 9];  col.funcfieldmaps=[10 11 12];  col.struct=13;  col.func_nscans=[14 15 16];  % log.datalog contents 
    scan=load([where.where filesep '1 Preprocessing' fs 'i_scanningdetails.mat']);
    
    % Log
    errorlog=cell(1,1); e=1;
    % Interface
    disp('===================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp('Requested analysis:'); disp(request); disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.subjects); end
    disp(' '); disp(['Data location: ' where.data_brain]); disp(' ')
    input('Hit Enter to start      ')
    disp('====================================')
    
end

%% Sort + Import 

if request.sort==1
    spm_jobman('initcfg');
    request.func=1;
    request.fieldmaps=1;
    request.structs=1; 
     
   
    for s=1:log.n_subjs
        % no need for any untarring. archive code for untarring: import_archive([ws.from spm_select('list', ws.from, [log.datalog{s+1,4} '.' num2str(log.datalog{s+1,4+i}) '.tar$'])],ws.to); ws.rename='java.io.file([ws.to log.datalog{s+1,4} ''.'' num2str(log.datalog{s+1,4+i})]).renameto(java.io.file([ws.to ''func_r'' num2str(i)]));'; eval(ws.rename)
        %         try
        disp(['Subject ' num2str(s) ' (' log.datalog{s+1} ')  --------- '])
        ws.from=[where.data_brain filesep log.datalog{s+1} filesep '0 Archive' filesep];
        ws.to=[where.data_brain filesep log.datalog{s+1} filesep '1 Preprocessed' filesep];
        if isdir(ws.to)==0; mkdir(ws.to); end
        
        % Compile index of scans: file name, session_no.-img_no., session_no
 
        ws.funcs= cell2mat(log.datalog(s+1,  col.funcsess));
        ws.scans=cellstr(spm_select('list', ws.from, '.*')); ws.scans= ws.scans(cellfun(@(x)~strcmp(x(1),'.'),  ws.scans(:,1)),:);  % This last bit is for mac running 
        ws.scans=[ws.scans   cellfun(@(x)x(strfind(x, 'PHOEBE')+7:strfind(x, '.2014.')-1), ws.scans,'uniformoutput',0)];
        ws.scans(:,3)= cellfun(@(x)str2num( x(1:strfind(x,'.')-1)),  ws.scans(:,2),'uniformoutput',0); ws.scans=sortrows(ws.scans,[3 2]);  ws.sc=cell2mat(ws.scans(:,3));
        scanlog=[{'Filename' 'Sess-Img' 'Session'}; ws.scans];  save([ws.to 'scanlog.mat'], 'scanlog'); scanlog=[];
        
        if request.func
            disp('1. Sorting Fxnals --- ')  
            for i=1:3
                ws.sn= ws.scans( ws.sc== log.datalog{s+1,  col.funcsess(i)});  
                matlabbatch{i}.spm.util.dicom.data=cellstr([repmat(ws.from, size(ws.sn,1),1) char(ws.sn)]);
                matlabbatch{i}.spm.util.dicom.root = 'flat';
                matlabbatch{i}.spm.util.dicom.outdir = {[ws.to 'Func_s' num2str(i)]}; mkdir([ws.to 'Func_s' num2str(i)])
                matlabbatch{i}.spm.util.dicom.convopts.format = 'img';
                matlabbatch{i}.spm.util.dicom.convopts.icedims = 0;
            end
            spm_jobman('run',matlabbatch); matlabbatch=[];
        end
        
        if  request.fieldmaps
            disp('2. Fieldmaps --- ')
            if strcmp(log.subjects{s}, 'p12_AL')==0
                ws.sn= ws.scans(find( (ws.sc== log.datalog{s+1,  col.funcfieldmaps(1)})+(ws.sc== log.datalog{s+1,  col.funcfieldmaps(1)}+1 ) ),1);
                matlabbatch{1}.spm.util.dicom.data=cellstr([repmat(ws.from, size(ws.sn,1),1) char(ws.sn)]);
                matlabbatch{1}.spm.util.dicom.root = 'flat';
                matlabbatch{1}.spm.util.dicom.outdir = {[ws.to 'Fieldmap']}; mkdir([ws.to 'Fieldmap'])
                matlabbatch{1}.spm.util.dicom.convopts.format = 'img';
                matlabbatch{1}.spm.util.dicom.convopts.icedims = 0;
                spm_jobman('run',matlabbatch); matlabbatch=[];
            end
        end
        
        if request.structs
            disp('3. Structurals --- ')
            ws.sn= ws.scans(find((ws.sc== log.datalog{s+1,  col.struct})),1);
            matlabbatch{1}.spm.util.dicom.data=cellstr([repmat(ws.from, size(ws.sn,1),1) char(ws.sn)]);
            matlabbatch{1}.spm.util.dicom.root = 'flat';
            matlabbatch{1}.spm.util.dicom.outdir = {[ws.to 'Structural']}; mkdir([ws.to 'Structural'])
            matlabbatch{1}.spm.util.dicom.convopts.format = 'img';
            matlabbatch{1}.spm.util.dicom.convopts.icedims = 0;
            spm_jobman('run',matlabbatch); matlabbatch=[];
        end
        
        %
        ws=[];
        %         catch
        %             errorlog{e,1}=['failed: sort & import --- ' log.datalog{s+1,1}];  disp(errorlog{e,1}); e=e+1;
        %         end
    end
end

% Delete dummy scans (func)
if request.deletedummyvolumes==1
    disp('############### step 2: deleting dummies ###############')
    for s=1:log.n_subjs
        disp(['Subject ' num2str(s) ' (' log.subjects{s} ')  --------- '])
        ws.where=[where.data_brain filesep log.datalog{s+1} filesep '1 Preprocessed' filesep];
        for sess=1:3
            ws.b.where=[ws.where 'Func_s' num2str(sess) fs'];
            for j=1:scan.nDummyVols
                f=spm_select('list', ws.b.where, ['^f.*00000' num2str(j) '-01']); 
                if isempty(f)==1; erorrlog{e,1}=['error. could not find dummy to delete  - ' log.subjects{s} ' -  session ' num2str(sess) '  dummy # ' num2str(j)];  end; e=e+1;
                for i=1:size(f,1), delete([ws.b.where f(i,:)]); end
            end
        end
            ws=[]; 
    end
end

    
%% Bias correction

if request.biascorrect==1,
    disp('########### bias correction ###########')
    spm_jobman('initcfg')
    for s=1:log.n_subjs
        disp(['Subject ' num2str(s) ' (' log.datalog{s+1} ')  --------- '])
        ws.to=[where.data_brain filesep log.datalog{s+1} filesep '1 Preprocessed' filesep];
        if  strcmp(log.datalog{s+1}, 'p12_AL')==1; ws.nb=2; else ws.nb=3;  end
        for i=1:ws.nb, 
            spm_biascorrect([ws.to filesep 'Func_s' num2str(i)]); 
        end
    end

end
             
 
%% END

disp('===================================='), w.c=clock;
disp(['END Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]), disp(' ')
disp('Analysis completed:'); disp(request)
disp(['No. of subjects: ' num2str(log.n_subjs)]); disp(' ')
disp('Error log?'); disp(errorlog); disp(' ')
disp('====================================')
try % Notify researcher
    f_sendemail('kurzlich', strcat('Analysis batchscript is complete (s1_SortImportBiascorrect)'), ' ',1);
end
