% Delete specific files/folders for all subjects 
clear all;close all hidden; clc

% Request specific
log.specificsubjects={}; % BLANK to process all subjects 

for o1=1:1 % General settings and specifications
    w=pwd;  if strcmp(w(1), '/')==1;  where.where='/Users/EleanorL/Dropbox/SCRIPPS/3 Pleasure purpose/';  where.data_brain='/Users/EleanorL/Desktop/3 PLPR';  
    else where.where='C:\Users\e.loh\Dropbox\SCRIPPS\3 Pleasure purpose\';  where.data_brain= 'D:\1 PLPP\1 MRI data'; 
    end  
    addpath(where.where), where.beh=[where.where filesep '3 Behaviour' filesep]; 
    if ischar(log.specificsubjects)==1 && strcmp(log.specificsubjects, 'AllDataOK'), log.specificsubjects= {'p01_YH';'p02_MI';'p03_AY';'p04_AW';'p05_CA';'p06_BT';'p07_HC';'p08_KC';'p09_KJ';'p10_YC';'p11_BI';'p13_MS';'p14_SK';'p15_MM';'p16_SH';'p17_BB';'p18_WB';'p19_HB';'p20_LZ';}; end ;
    [n t r]=xlsread([where.beh 'datalog_plpr.xlsx']); 
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(r, log.specificsubjects,  r, 'All'); 
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
    disp('Requested analysis: ADHOC SCRIPT'); disp(' ')
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0;  if sum(strcmp(log.subjects, {'p01_YH';'p02_MI';'p03_AY';'p04_AW';'p05_CA';'p06_BT';'p07_HC';'p08_KC';'p09_KJ';'p10_YC';'p11_BI';'p13_MS';'p14_SK';'p15_MM';'p16_SH';'p17_BB';'p18_WB';'p19_HB';'p20_LZ';}) ) ==19;  disp('   Subjects w complete data only (n=19)');  else disp('   Subset of subjects only:'); disp(log.specificsubjects); end;  end; disp(' ')
    disp(' '); input('Hit Enter to start      ')
    disp('=======================================================')
    
end

mods= {
    'c2_ChoUnchoPLnPP'
    'c2b_ChoUnchoPLnPP'
    'c3_PLnPP_conflict'
    'c3b_PLnPP_conflict'
    'cl2b_ChoUnchoPLnPP'
    'cl3b_PLnPP_conflict'
};
%%

for s=1: log.n_subjs
    ws.c=clock;  disp(['Subject ' num2str(s) '   -  ' log.subjects{s} '   [' num2str(ws.c(4)) ':' num2str(ws.c(5)) ']  ------------------']);
    ws.subfol=[where.data_brain filesep log.subjects{s} filesep]; 
    ws.subfol=[ws.subfol '2 First level' filesep];
%%
 
cd(ws.subfol)
for m=1:length(mods)
try 
%     rmdir([mods{m} ' Contrasted'], 's')
    
    delete([log.subjects{s} '_onsets_' mods{m} '.mat'])
        
    
catch; disp(['dd ' num2str(m)])
end   
end
%%
%         eval('java.io.File(ws.old).renameTo(java.io.File(ws.new));')
ws=[];
end

% openvar subids

error('DONE!')

try % Notify researcher
    f_sendemail('kurzlich', strcat('DONE with ad hoc script'), ' ',1);
end


