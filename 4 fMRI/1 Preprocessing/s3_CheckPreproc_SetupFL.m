% Plot movement, Organize, Checkreg
% clear all; close all hidden; clc
clear all; clc

% Requested analysis
request.plotmovement=0;
request.identifymovescans=0;
request.checkreg=1;
request.setup_FLfol=0;
request.preproc_prefix='s8wub'; 
log.specificsubjects={
%     'p01_YH';'p02_MI';'p03_AY';'p04_AW';'p05_CA';'p06_BT';'p07_HC';'p08_KC';'p09_KJ';'p10_YC';
%     'p11_BI';'p12_AL';'p13_MS';'p14_SK';'p15_MM';'p16_SH';'p17_BB';'p18_WB';'p19_HB';'p20_LZ'
    };
for o1=1:1 % General settings and specifications    
    
    % Load subjects
    w=pwd;  if strcmp(w(1), '/')==1;  where.where='/Users/EleanorL/Dropbox/SCRIPPS/3b PLPP fmri';  where.data_brain='/Users/EleanorL/Desktop/3 PLPR/1 Brain data';  
    else where.where='C:\Users\e.loh\Dropbox\SCRIPPS\3b PLPP fmri';  where.data_brain= 'D:\1 PLPP\1 MRI data';  where.spm='C:\toolbox\64\spm8';
    end   
    [n t log.datalog]=xlsread([where.data_brain fs  'datalog_plpr.xlsx']); path(pathdef), addpath(where.where)
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    col.funcsess=[7 8 9];  col.funcfieldmaps=[10 11 12];  col.struct=13;  col.func_nscans=[14 15 16];  % log.datalog contents 
    scan=load([where.where filesep '1 Preprocessing' fs 'i_scanningdetails.mat']);
    errorlog=cell(1,1); e=1;  % Log 
    f_mat=@(A,x)A(x); 
    
    % Interface
    disp('===================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ])
    disp(' '); disp('Requested:'); disp(request)
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.specificsubjects); end
    disp(' '); disp(['Data location: ' where.data_brain]); disp(' ')
    input('Hit Enter to start      ')
    disp('====================================')
    
end


%% (1) Plot movement regressors
 
if request.plotmovement==1 || request.identifymovescans
    disp('----------------- Collecting movement parameters ----------------- ')
    movement=cell(log.n_subjs,7); checkscan =  repmat( {[]}, log.n_subjs,3); movediff= movement;
    request.badmove=[1.5  1];  % mm translation, degree rotations 
    for s=1:log.n_subjs
        disp(['Subject ' num2str(s) '  (' log.subjects{s} ')']) 
        for b=1:3
            f=spm_select('List', [where.data_brain filesep log.subjects{s} filesep  '1 Preprocessed' filesep 'Func_s' num2str(b)], '.txt$'); 
            if isempty(f)==1, movement{s,b}=zeros(10,6) ; disp(['     Where''s your movement reg file? b' num2str(b)]);  
            else   movement{s,b}=load([where.data_brain filesep log.subjects{s} filesep  '1 Preprocessed' filesep 'Func_s' num2str(b) filesep f(1,:)]);  
            end
            
            % Movement = difference from previous scan 
            movement{s,b}(:, 4:6)= movement{s,b}(:, 4:6)*180/pi;  % Convert to degrees from radians 
            movediff{s,b}= movement{s,b} -  [movement{s,b}(1, :); movement{s,b}(1:end-1, :) ]; 
             
             % Find badduns
             checkscan{s,b}= scan.nDummyVols + find(sum(movediff{s,b} >[repmat(request.badmove(1),size(movediff{s,b},1),3) repmat(request.badmove(2),size(movediff{s,b},1),3)],2)>0);   
             disp([     'B' num2str(b) ': ' num2str( length(checkscan{s,b})) ' to check'])
        end
        movement{s,4}=log.subjects{s};
    end
    
    if request.plotmovement==1 % Plot
        f.plotcols=3;  f.figwidth= 1300; f.figheight=1400; f.fontsize=15; f.fontsize_title=30;
        f.markesize= 4; 
        f.subplot_VerHorz=[0.03 0.03]; f.fig_BotTop=[0.05 0.05]; f.fig_LeftRight=[0.1 0.03]; k=1;
        figure('Name', 'Movement parameters', 'Position', [400 100 1200 1000], 'color', 'w')
        for s=1:log.n_subjs
            for i=1:3
                subtightplot( round(log.n_subjs*3/f.plotcols),  f.plotcols, k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight); k=k+1; 
                scatter(1:size(movediff{s,i},1),  movediff{s,i}(:,1), f.markesize, [0.8 0 0]), hold on
                scatter(1:size(movediff{s,i},1),  movediff{s,i}(:,2), f.markesize, [0.8 0 0]), hold on
                scatter(1:size(movediff{s,i},1),  movediff{s,i}(:,3), f.markesize, [0.8 0 0]), hold on
                scatter(1:size(movediff{s,i},1),  movediff{s,i}(:,4), f.markesize, [0 0 0.8]), hold on
                scatter(1:size(movediff{s,i},1),  movediff{s,i}(:,5), f.markesize, [0 0 0.8]), hold on
                scatter(1:size(movediff{s,i},1),  movediff{s,i}(:,6), f.markesize, [0 0 0.8]), hold on 
                plot(1:size(movediff{s,i},1), ones(size(movediff{s,i},1),1)*request.badmove(1),'k');  plot(1:size(movediff{s,i},1), ones(size(movediff{s,i},1),1)*-request.badmove(1),'k')
                plot(1:size(movediff{s,i},1), ones(size(movediff{s,i},1),1)*request.badmove(2),'k');  plot(1:size(movediff{s,i},1), ones(size(movediff{s,i},1),1)*-request.badmove(2),'k')
                xlim([1 size(movediff{s,i},1)])
                if s==1, title(['s' num2str(i)],'FontSize',   f.fontsize_title ), end
                if i==1, ylabel(log.subjects{s},'FontSize',   f.fontsize), end
                ylim([-request.badmove(1)*1.5 request.badmove(1)*1.5])
            end
        end
    end
end
 

%% (2) CheckReg 

if request.checkreg==1 
    disp('--------------- Collecting scans for CheckReg --------------------')
    w.nScansPerFigure=8;
    %
    images=cell(log.n_subjs*2,1); i=1;
    for s=1:log.n_subjs
        ws.preprox_prefix =  request.preproc_prefix; 
        if  strcmp(log.subjects{s}, 'p12_AL'); ws.preprox_prefix=strrep(ws.preprox_prefix, 'u', 'r');  ws.nb=2; 
        else  ws.nb=3; 
        end 
        
        for b=1:ws.nb 
            wb.where=[where.data_brain filesep log.subjects{s} filesep '1 Preprocessed' filesep 'Func_s' num2str(b) filesep];
            
            % Identify scans
            f=spm_select('List', wb.where, ['^' ws.preprox_prefix  '.*img$']);
            if  isempty(f); disp(['Error: no scans for ' log.subjects{s} ' block ' num2str(b)]); end
            wb.runids = unique(cellstr(f(:, 1: f_mat(strfind(f(1,:), '-'),2)-1)));    % 1 scan per run
            images{i,1}=[wb.where char(f_mat(cellstr(spm_select('List', wb.where, ['^' wb.runids{1}(1,:) '.*img$'])), randi(size(spm_select('List', wb.where, ['^' wb.runids{1}(1,:) '.*img$']),1)))) ',1'];
            images{i,2}=[log.subjects{s} ' s' num2str(b) ];  i=i+1; 
        end
    end
    
    % Collate into figures
    f=1; i=1; figs=cell(ceil(size(images,1)/w.nScansPerFigure),1);
    for d=1:size(images,1)
        figs{f}{i,1}=char(images{d,1});
        images{d,3}=f;
        if i==w.nScansPerFigure i=1; f=f+1;
        else i=i+1; 
        end
    end
    
    % Display / Instructions for display
    if w.nScansPerFigure>=size(images,1)
        matlabbatch{1}.spm.util.checkreg.data=figs{1};
        spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);
    else
        disp(' --------------------- INSTRUCTIONS FOR CHECKREG -------------------------------------')
        disp(['Multiple figures to display (' num2str(size(figs,1)) ' figures)'])
        disp('Note: To change no. of images per display, specify in script (at beginning of module)')
        disp('See variable ''images'' to identify subject scans (Col 3)'), disp(' ')
        disp('To display each batch, specify value of ''f'' and execute following command:'), disp(' ')
        disp('f = # ; eval(checkregcommand) '), disp(' ')        
        disp(' -----------------------------------------------------------------------------------------------')
        spm_jobman('initcfg');  checkregcommand='matlabbatch{1}.spm.util.checkreg.data=figs{f}; spm_jobman(''run'' , matlabbatch);';
    end
end

 



%% (3) Prep for First level

if request.setup_FLfol==1
    
    request.transfer_funcs=1;
    disp('--------------- Preparing folder for 1st level --------------------')
    for s=1:log.n_subjs
        disp(['Subject  ' num2str(s) '  (' log.subjects{s} ') '])
        wb.where=[where.data_brain filesep log.subjects{s} filesep];
        wb.wherefrom=[wb.where '1 Preprocessed' filesep];
        wb.whereto=[wb.where '2 First level' filesep];
        if isdir(wb.whereto)==0; mkdir(wb.whereto); mkdir([wb.whereto 'Preproc func s1']),  mkdir([wb.whereto 'Preproc func s2']),  mkdir([wb.whereto 'Preproc func s3']), end 
        if  strcmp(log.subjects{s}, 'p12_AL'),    ws.preprox_prefix =[strrep(request.preproc_prefix , 'u', []) '*'];ws.nb= 2;  else  ws.preprox_prefix =request.preproc_prefix ; ws.nb= 3;  end

        if request.transfer_funcs
            for r=1:ws.nb % Collect functionals
                cd([wb.wherefrom 'Func_s' num2str(r)])
                
                wb.files=spm_select('List', pwd, ['^' ws.preprox_prefix  '.*']);
                disp([log.subjects{s} '  run ' num2str(r) ' : '  num2str(length(wb.files)) ' files'])
                for i=1:size(wb.files,1)  % Copy func files if not done yet
                    copyfile([wb.wherefrom 'Func_s' num2str(r) filesep wb.files(i,:)],[wb.whereto 'Preproc func s' num2str(r) filesep wb.files(i,:) ]);
                end
%                 disp('Transfer of func files turned off!')
            end
        end
        
        % Movement regressors   
        f=spm_select('List', [wb.wherefrom 'Func_s1'], '^rp.*txt$');    % Rating session
        if isempty(f);  disp(['    ' log.subjects{s} ' 1rating movement regs: failed'])
        else  copyfile( [wb.wherefrom 'Func_s1' filesep f], [wb.whereto log.subjects{s} '_1rating_regmovement.txt'])
        end
        f=spm_select('List', [wb.wherefrom 'Func_s2'], '^rp.*txt$');     % Choice main session 
        if isempty(f);  disp(['    ' log.subjects{s} ' 2choice movement regs: failed'])
        else  copyfile( [wb.wherefrom 'Func_s2' filesep f], [wb.whereto log.subjects{s} '_2choice_regmovement.txt'])
        end 
        if  ~strcmp(log.subjects{s}, 'p12_AL'),     
            f=spm_select('List', [wb.wherefrom 'Func_s3'], '^rp.*txt$');     % Choice Lab
            if isempty(f);  disp(['    ' log.subjects{s} ' 3choice lab movement regs: failed'])
            else  copyfile( [wb.wherefrom 'Func_s3' filesep  f], [wb.whereto log.subjects{s} '_3choicelab_regmovement.txt'])
            end
            
            % Choice: concatenate main and in-lab sessions
            wb.rmotion = [];  f=spm_select('List', [wb.wherefrom 'Func_s2'], '^rp.*txt$');
            wb.b=load([wb.wherefrom 'Func_s2' filesep f(1,:)]);  
            wb.rmotion=  [wb.rmotion; [wb.b wb.b(:,1)*0+1]];
            f=spm_select('List', [wb.wherefrom 'Func_s3'], '^rp.*txt$');
            wb.b=load([wb.wherefrom 'Func_s3' filesep f(1,:)]);  
            wb.rmotion=  [wb.rmotion; [wb.b abs(wb.b(:,1)*0)]]; 
            [wb.printok]=print2txt(wb.whereto, [log.subjects{s} '_4choiceall_regmovement'], wb.rmotion);
        end
        
        wb=[]; ws=[]; 
    end
end

%%  

disp(' ###################################################'), disp(' DONE ' )
if request.checkreg==1, disp('CHECKREG: See command window (earlier command) for display instructions'), end
disp(' ###################################################')



