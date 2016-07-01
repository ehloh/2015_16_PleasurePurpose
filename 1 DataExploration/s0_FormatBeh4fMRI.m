% s0_FormatBeh4fMRI
clear all, close all, clc 



req.save = 0; 

for o=1:1
    where.where = '/Users/EleanorL/Dropbox/SCRIPPS/3 Pleasure purpose'; 
%     where.where ='C:\Users\e.loh\Dropbox\SCRIPPS\3 Pleasure purpose';
    where.data = [where.where fs '3 Behaviour' fs ]; 
    cd(where.data), subs= dir('p*'); subs=cellstr(char(subs.name)); 
    path(pathdef);addpath([where.where filesep '4a Beh analysis basic'])
end 

%% Definitions 

% Binning
log.define.scorebins={[1 2 3]; [4 5 6 7];  [8 9 10]};

% Conflict formulation
%   Raw conflict scores are always calculated in main
%   script, whereas binned scores are calc'd in functions 
log.define.conflict ='Max - abs diff';
log.fcf_raw= @(x,y)((2*10-1)- abs(x- y) );
log.fcf_binned= @(x,y)(   ( 2*size(log.define.scorebins,1)-1 )    - abs(x- y) );
binfxns.fcf = log.fcf_binned;   % Conflict formulation for binned scores 

%%

for o=1:1  % Original columns (read out from excel sheet)
    % Session 1
    col1.Trialnum=1;
    col1.Eventnum=2;
    col1.PL=3;
    col1.PP=4;
    col1.Onset_event=5;
    col1.Duration_event=6;
    col1.Onset_PLrate=7;
    col1.RT_PL=8;
    col1.Onset_PPrate=9;
    col1.RT_PP=10;
    col1.Onset_PLmotor=11;
    col1.Onset_PPmotor=12;
    col1.Duration_PLmotor=13;
    col1.Duration_PPmotor=14;
    col1.rateVivid=15;
    col1.rateValence=16;
    col1.rateArousal=17;
    col1.rateExperience=18;
    col1.rateHappy=19;
    col1.rateSocial=20;
    
    % Session 2
    col2.Trialnum=1;
    col2.Option1=2;
    col2.Option2=3;
    col2.Onset_Options=4;
    col2.RT_Choice=5;
    col2.Onset_Choice=6;
    col2.Duration_Choice=7;
    col2.Choice=8;
    col2.PL1=9;
    col2.PL2=10;
    col2.PP1=11;
    col2.PP2=12;
    
    % Session 3
    col3.Trialnum=1;
    col3.Option1=2;
    col3.Option2=3;
    col3.Onset_Options=4;
    col3.RT_Choice=5;
    col3.Onset_Choice=6;
    col3.Duration_Choice=7;
    col3.Choice=8;
    col3.PL1=9;
    col3.PL2=10;
    col3.PP1=11;
    col3.PP2=12;
end
for o=1:1  % New columns 
      
    req.col_rate = {
        'PLPP';  'TrialOK';
        'PLbin'; 'PPbin';'PLPPbin';  
        }; 
    req.col_choose= {
        'PLPP1'; 'PLPP2'; 'PLcf';'PPcf';'PLPPcf'; 'Session'; 
         'choPL';'choPP';'choPLPP';'unchoPL';'unchoPP';'unchoPLPP';
         'conchoPL';'conchoPP';'conchoPLPP';'Trialcf';'PPbest';'PLbest';'PLPPbest'; 
        'PLvcf';'PPvcf';'PLPPvcf';
        'd1m2.PL';'d1m2.PP';'d1m2.PLPP'; 'ad1m2.PL';'ad1m2.PP';'ad1m2.PLPP'; 
        'PLdraw';'PPdraw';'PLPPdraw'; 'chobest'; 'TrialOK'; 
        'marPL';'marPP';'marPLPP';'pChoPL';'pChoPP';'pChoPLPP';
        'RTchoPL';'RTchoPP';'RTchoPLPP';
        'Duration_Options'; 'Duration_OptionsChoice'; 'Duration_Trial';
        'PLbinvcf';'PPbinvcf';'PLPPbinvcf'; 'PLbin1';'PLbin2';'PPbin1';'PPbin2';'PLPPbin1';'PLPPbin2';'PLbincf';'PPbincf';'PLPPbincf'; 
        'd1m2.PLbin';'d1m2.PPbin';'d1m2.PLPPbin'; 'ad1m2.PLbin';'ad1m2.PPbin';'ad1m2.PLPPbin';
        'choPLbin';'choPPbin';'choPLPPbin';'unchoPLbin';'unchoPPbin';'unchoPLPPbin';
        'PLmarbin'; 'PPmarbin';'PLPPmarbin'; 'pChoPLbin';'pChoPPbin';'pChoPLPPbin';
    }; 
    col2.PL=[col2.PL1 col2.PL2]; col3.PL=[col3.PL1 col3.PL2]; col2.PP=[col2.PP1 col2.PP2]; col3.PP=[col3.PP1 col3.PP2];
    for i=1:length(req.col_rate),  eval(['col1.'  req.col_rate{i} '=structmax(col1)+1;']), end
    for i=1:length(req.col_choose),  eval(['col2.'  req.col_choose{i} '=structmax(col2)+1;']), eval(['col3.'  req.col_choose{i} '=structmax(col3)+1;']), end;   
    col2.PLPP=[col2.PLPP1 col2.PLPP2]; col3.PLPP=[col3.PLPP1 col3.PLPP2];
    col1.PLrating = col1.PL;  col1.PPrating = col1.PP;  col2.PLrating = col2.PL;  col2.PPrating = col2.PP;  col3.PLrating = col3.PL;  col3.PPrating = col3.PP; 
    
end 



for s=1:length(subs)
    cd([where.data subs{s}])  
    disp(subs{s})
    
    %% Session 1 
    [n t r] = xlsread([subs{s} '_behdata.xlsx'], 'session1_all');  
    data= n(1:80, 1:20); 
    data(:, col1.PLPP)= data(:, col1.PL) + data(:, col1.PP);
    data(:, col1.TrialOK)=1; 
    data(isnan(data(:, col1.PP)) | isnan(data(:, col1.PP)), col1.TrialOK)=0;  
    for i=1:length(log.define.scorebins)
        for ii=1:length(log.define.scorebins{i}) 
            data(data(:, col1.PL)== log.define.scorebins{i}(ii), col1.PLbin)=i;
            data(data(:, col1.PP)== log.define.scorebins{i}(ii), col1.PPbin)=i; 
        end 
    end 
    data(:, col1.PLPPbin)= data(:, col1.PLbin) + data(:, col1.PPbin);
    %
    col=col1; if req.save , save([subs{s} '_1Rating.mat'], 'data','col', 'log'), end 
    data1=data; col=[];
      
    %% Session 2
    [n t r] = xlsread([subs{s} '_behdata.xlsx'], 'session2_all');  
    data= n(1:80, 1:12); 
    data(data(:, col2.Choice)==-1, col2.Choice)=nan;   
    data( :, col2.Choice)=data( :, col2.Choice)+1;   %  In excel, 0= Choice 1, 1= Choice 2 
    data = f_calcfromraw(data, col2); 
    data(:, col2.PLcf)= log.fcf_raw(data(:, col2.PL1),  data(:, col2.PL2));
    data(:, col2.PPcf)= log.fcf_raw(data(:, col2.PP1),  data(:, col2.PP2));
    data(:, col2.PLPPcf)= log.fcf_raw(data(:, col2.PLPP1),  data(:, col2.PLPP2));   
    %
    data = f_calcnew(data, col2,log.define.scorebins, binfxns); 
    %
    col=col2;  if req.save, save([subs{s} '_2Choice.mat'], 'data','col', 'log'), end
    data2=data; col=[];
    
    %% Session 3
    if strcmp(subs{s}, 'p12_AL')==0 
        [n t r] = xlsread([subs{s} '_behdata.xlsx'], 'session3_all');  
        data= n(1:20, 1:12); 
        data(data(:, col3.Choice)==-1, col3.Choice)=nan; 
        data( :, col3.Choice)=data( :, col3.Choice)+1;   %  In excel, 0= Choice 1, 1= Choice 2  
        data = f_calcfromraw(data, col3);        
        data(:, col3.PLcf)= log.fcf_raw(data(:, col3.PL1),  data(:, col3.PL2));
        data(:, col3.PPcf)= log.fcf_raw(data(:, col3.PP1),  data(:, col3.PP2));
        data(:, col3.PLPPcf)= log.fcf_raw(data(:, col3.PLPP1),  data(:, col3.PLPP2));
        %
        data = f_calcnew(data, col3,log.define.scorebins, binfxns); 
        %
        col=col3; if req.save, save([subs{s} '_3ChoiceLab.mat'], 'data','col', 'log'), end
        data3=data; col=[];
    else  data3=[];  
    end 
    
    for o=1:1    % Manually check that values are ok?
        docheck = 1; 
        if docheck 
%             d2=data, col2=col;
%             d2(d2(:, col2.Choice)==1, :)
            d2= data2;  
%             d2(d2(:, col2.Choice)==1, [col2.choPL col2.choPP col2.unchoPL col2.unchoPP])-d2(d2(:, col2.Choice)==1, [col2.PL1 col2.PP1 col2.PL2 col2.PP2 ])  % Session 2, Chosen & unchosen logged properly?
            
            
%             
%             if s==2; er; end 
%             
%             d3=data3; 
%              
%             
%             d3(:, col3.choPLbin)-d3(:, col3.choPPbin)
%             
%             openvar col3
%             d3
%             
%             
            
            
        end
    end
    
    
    if ~strcmp(subs{s}, 'p12_AL')
        sdd(1:size(data3,1), s)=  data3(:, col3.choPLbin)- data3(:, col3.choPPbin);
    end
     
    % Save 
    if req.save, save([subs{s} '_behdata.mat'], 'data1','col1', 'data2','col2', 'data3','col3', 'log'), end
end 
if req.save, disp('SAVED'); else disp('NO SAVING!'); end 

openvar sdd

 