% Look at regressors and check for correlatedness
clear all;close all hidden; clc, f_mat=@(A,x)A(x);

% Requested analysis
log.specificsubjects={
%     'p01_YH'; 'p02_MI';
%     'p03_AY';'p04_AW';'p05_CA';'p06_BT';'p07_HC';'p08_KC';'p09_KJ';'p10_YC';'p11_BI';
% 'p12_AL'; 'p13_MS';'p14_SK';'p15_MM';'p16_SH';'p17_BB';'p18_WB';'p19_HB';'p20_LZ' 
    };


% Which model? ########################
% log.onsetsmodel='Beh2';  
% log.onsetsmodel='r1_PLnPPrating'; 
log.onsetsmodel='c2_ChoUnchoPLnPP';  
% log.onsetsmodel='c3b_PLnPP_conflict';  

 

% SUBSAMPLE ?? 
log.specificsubjects= f_subsample('beh1','PLvPP_abs_r04'); 
% log.specificsubjects= f_subsample(log.onsetsmodel,'regr_ot_PLvPP_abs_r04'); 
 


for o1=1:1 % General settings and specifications 
    
    % Load subjects
    w.w=pwd; 
    if strcmp(w.w(1), '/')==1;  where.where='/Users/EleanorL/Dropbox/SCRIPPS/3b PLPP fmri';  
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
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    if  strcmp(log.onsetsmodel(1:2), 'cl') + strcmp(log.onsetsmodel, 'Beh3')>0
        log.specificsubjects= log.subjects(~strcmp(log.subjects, 'p12_AL'));  disp('Excluding p12_AL from this model');  
        [log.subjects log.n_subjs log.datalog] = f_selectsubjects(log.datalog, log.specificsubjects, [log.datalog vertcat('include_all', num2cell(ones(size(log.datalog,1)-1,1)))], 'include_all');
    end  
    
    
    % Stable settings  
    log.FirstLevelThread=[]; log.prefix='s8wubf';  
    f_mat=@(A,x)A(x); errorlog=cell(0,1); e=1; 
    
    % What sort of model is this? 
    switch log.onsetsmodel(1: regexp(log.onsetsmodel, '\d')-1)
        case 'r',  log.modeltype= 'Rating';  
        case 'c',  log.modeltype= 'Choice';  
        case 'cl',  log.modeltype= 'ChoiceLab';   
    end   
 

    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0;disp('   Subset of subjects only:'); disp(log.subjects);  end
    disp(' '); disp(['Data location (brain): ' where.data_brain])
    disp(' ');disp('CHECK HERE: 1st and 2nd level models ---------'); disp(' ')
    disp(['             First level model:  ' log.onsetsmodel])
    input('Hit Enter to start      ')
    disp('=======================================================')
    
end
where.data_brain = '/Users/EleanorL/Dropbox/WorkPC/PLPP'; 
where.data_beh = '/Users/EleanorL/Dropbox/SCRIPPS/3 Pleasure purpose/3 Behaviour';

%% Which Variables are included in this model? (sample 1st subject)
 
if strcmp(log.onsetsmodel(1:3), 'Beh') 
    log.session = log.onsetsmodel(4);  
    switch log.session
        case '1', log.VarName={'PL';'PP'};   
        case '2', log.VarName={'PL1'; 'PL2'; 'PP1'; 'PP2'; 'PLcf'; 'PPcf'; 'choPL'; 'choPP'; 'unchoPL'; 'unchoPP';};   
            
            log.VarName={'choPL'; 'choPP'; 'unchoPL'; 'unchoPP';};   
            
        case '3', log.VarName={'PL1'; 'PL2'; 'PP1'; 'PP2'; 'PLcf'; 'PPcf'; 'choPL'; 'choPP'; 'unchoPL'; 'unchoPP';};   
        otherwise, error('Session invalid');
    end 
    colcommand =   strjoin(cellfun(@(x)['col.' x ' '], log.VarName, 'uniformoutput',0));
    d_var= cell(log.n_subjs,1);
    
    for s=1:log.n_subjs % counting scans
        disp(['Subject ' num2str(s) '   -  ' log.subjects{s}])
        
        % Load + sort behavioural data
        ws.f= load([where.data_beh fs log.subjects{s} fs log.subjects{s} '_behdata.mat']); 
        eval(['col=ws.f.col' num2str(log.session) ';'])
        eval(['ws.d =ws.f.data' num2str(log.session) ';']) 
        eval(['d_var{s} =ws.d(:,[' colcommand ']);'])  
    end 
    
    
else  % Read out conditions of interest in this fMRI model (Assuming all subject have the same of-interest regressors)
        %               log.FLregs: col 1=regressor name, col 2= regressor num 
        try wc=load([where.data_brain filesep log.subjects{1} filesep '2 First level' log.FirstLevelThread filesep log.onsetsmodel ' Contrasted' filesep 'SPM.mat']);
        catch; wc=load([where.data_brain filesep log.subjects{1} filesep '2 First level' log.FirstLevelThread filesep log.onsetsmodel ' Estimated' filesep 'SPM.mat']);
        end; SPM = wc.SPM; 
        log.FLregs =  f_mat(   cellfun(@(x)x(7:length(x)), wc.SPM.xX.name(1:3:end)','uniformoutput',0), 1:find(cellfun(@(x)strcmp(x(1:2), 'n_'), cellfun(@(x)x(7:length(x)), wc.SPM.xX.name(1:3:end)','uniformoutput',0)),1,'first')-1);   % Up to b4 1st no-interest regressor
        log.FLregs(:,2) =    cellfun(@(x, m)find(strcmp(x,m)), log.FLregs, repmat({cellfun(@(x)x(7:length(x)),  wc.SPM.xX.name ,'uniformoutput',0)}, length(log.FLregs),1), 'uniformoutput',0); 
        disp('FL regressors (of-interest only):');  disp(log.FLregs) 
        
        % Get subject regs 
        d_var= cell(log.n_subjs,1); 
        for s=1:log.n_subjs  
            try wc=load([where.data_brain filesep log.subjects{s} filesep '2 First level' log.FirstLevelThread filesep log.onsetsmodel ' Contrasted' filesep 'SPM.mat']);
            catch; wc=load([where.data_brain filesep log.subjects{s} filesep '2 First level' log.FirstLevelThread filesep log.onsetsmodel ' Estimated' filesep 'SPM.mat']);
            end
            d_var{s,1} =  wc.SPM.xX.X(:,  cell2mat(log.FLregs(:,2) )); 
        end
         
        % Output  
        log.VarName=    strrep( strrep(log.FLregs(:,1),'*bf(1)','') ,'^1',''); 
end
log.nVar= length(log.VarName);

%% Correlation analysis 

% Execute correlations 
d_corr= nan(log.n_subjs,log.nVar^2); d_corrbeh=cell(log.nVar,log.nVar);
for s=1:log.n_subjs
%     d_corr(s,:) = f_mat(corr(d_var{s,1} ), 1:log.nVar^2);
    
    d_corr(s,:) = f_mat(corrcoef(d_var{s,1},'rows','pairwise'), 1:log.nVar^2);
    
    
end


% Req
req.r_limit= 0.4; 


close all hidden
f.plotcols=2;  f.figwidth= 800; f.figheight=800; f.fontsize=20; f.fontsize_title=30;
f.subplot_VerHorz=[0.05 0.10]; f.fig_BotTop=[0.05 0.05]; f.fig_LeftRight=[0.2 0.1];
f.markersize= 15; 
for o=1:1 % Plots 
    figure('Name', 'Pairwise correlatedness', 'NumberTitle', 'off', 'Position', [200 150 f.figwidth f.figheight], 'Color', 'w');  k=1;

    % Group plot
    subtightplot(log.nVar,  f.plotcols,((1:log.nVar)-1)*f.plotcols+1,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);
    %   subtightplot(log.nVar+1,  f.plotcols,k ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);  k=k+1;
    imagesc( reshape(mean(d_corr), [log.nVar log.nVar])), axis square,  colorbar
    set(gca,'xtick', 1:log.nVar,'ytick', 1:log.nVar,'xticklabel', 1:log.nVar, 'yticklabel',  cellfun(@(x,n)[x '  [' num2str(n) ']   '], log.VarName, num2cell(1:log.nVar)','uniformoutput',0), 'FontSize',f.fontsize)
    title('Mean correlation r', 'FontSize', f.fontsize_title )
    
    for v=1:log.nVar   % Group plot w SEs
        subtightplot(log.nVar,  f.plotcols, (v-1)*f.plotcols+2 ,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);
        %   subtightplot(log.nVar+1,  f.plotcols, k,f.subplot_VerHorz,f.fig_BotTop, f.fig_LeftRight);  k=k+1;
        barwitherr(f_mat(reshape(std(d_corr), [log.nVar log.nVar])/sqrt(log.n_subjs),    (v-1)*log.nVar+(1:log.nVar)),         f_mat(reshape(mean(d_corr), [log.nVar log.nVar]),    (v-1)*log.nVar+(1:log.nVar))    ,'y'   ); hold on
        scatter(f_mat(  repmat(1:log.nVar, log.n_subjs,1),   1:log.nVar*log.n_subjs),  f_mat(  d_corr(:, (1:log.nVar)+(v-1)*log.nVar),    1:log.nVar*log.n_subjs), f.markersize), hold on
        plot(0:log.nVar+1, ones(log.nVar+2,1)*req.r_limit, 'g'), hold on
        plot(0:log.nVar+1, ones(log.nVar+2,1)*-req.r_limit, 'g')
        xlim([0 length(log.VarName)+1]) 
        ylabel(log.VarName{v}, 'FontSize', f.fontsize)
        set(gca,'FontSize',f.fontsize)
        
        % Collect data in a more identifiable context  
        wv.vd=  d_corr(:, (1:log.nVar)+(v-1)*log.nVar);
        for vv=1:log.nVar   
            d_corrbeh{v,vv}= wv.vd(:, vv);
        end 
        
    end
end

 
vrowcol=[2 3];
sd= sortrows([log.subjects num2cell(d_corrbeh{vrowcol(1),vrowcol(2)})],2);
sd_abs= sortrows([log.subjects num2cell(abs(d_corrbeh{vrowcol(1),vrowcol(2)}))],2);
openvar sd, openvar sd_abs
 [sd(:,1) sd_abs(:,1)]
 
 
 
%% 
 
 
 