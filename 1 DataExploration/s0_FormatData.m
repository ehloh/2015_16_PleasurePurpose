% Compile Phoebe's excel files into matlab data. Outputs: subjdata for behavioural modelling 
%
%   Note: This script allows for flexible changing of things (e.g. bins).
%           Files are thus always dated 
%   s0_FormatData4fMRI is for the final one (i.e. scores that are carried
%   thru to all fMRI analysis etc).
clear all;close all hidden; clc

request.do_mainstudy=1;
request.do_pilot=0; 

% [ Bin scores ? ] #########
log.scorebins=[];  request.filename_suffix=[]; % No binning 
% log.scorebins={[1 2 3 4]; [5  6 7];  [8 9 10]};   request.filename_suffix=' b3'; 
%     log.scorebins={[1 2]; [3 4]; [5 6]; [7 8]; [9 10]}; request.filename_suffix=' b5'; 

% log.scorebins={[1 2 3]; [4 5  6 7];  [8 9 10]};   request.filename_suffix=' b3'; 

% Request specific
log.specificsubjects={}; % BLANK to process all subjects
% log.specificsubjects='AllDataOK';  % Only complete data ok

for o1=1:1 % General settings and specifications
   
    % Load subjects
    where.where='/Users/EleanorL/Dropbox/SCRIPPS/3 Pleasure purpose';  
    where.beh=[where.where filesep '3 Behaviour']; 
    addpath(where.where)
    if ischar(log.specificsubjects)==1 && strcmp(log.specificsubjects, 'AllDataOK'), log.specificsubjects= {'p01_YH';'p02_MI';'p03_AY';'p04_AW';'p05_CA';'p06_BT';'p07_HC';'p08_KC';'p09_KJ';'p10_YC';'p11_BI';'p13_MS';'p14_SK';'p15_MM';'p16_SH';'p17_BB';'p18_WB';'p19_HB';'p20_LZ';}; end ;
    [n t r]=xlsread([where.beh filesep 'datalog_plpr.xlsx']); 
    [log.subjects log.n_subjs log.datalog] = f_selectsubjects(r, log.specificsubjects,r, 'All'); 
    if isempty(log.scorebins); log.define={'Scorebins'  'Unbinned' ' '};  else log.define={'Scorebins'  size(log.scorebins,1) log.scorebins};  end
    
    % Interface
    disp('======================================================='); w.c=clock;
    disp(['START Timestamp: ' date ' ' num2str(w.c(4)) ':' num2str(w.c(5)) ' hrs' ]); disp(' ')
    disp(['No. of subjects: ' num2str(log.n_subjs)])
    if isempty(log.specificsubjects)==0; disp('   Subset of subjects only:'); disp(log.specificsubjects); end; disp(' ')
    if isempty(log.scorebins)==1; disp('No score binning requested'); else disp(['Requested binning:   ' num2str(length(log.scorebins)) ' ordinal bins']); end
    disp(' '); input('Hit Enter to start      ')
    disp('=======================================================')
    
end

%%

for o=1:1  % Formulations  
    
    
    % Conflict formulation
	%   note that PLPP and PL/PP cf scores are equivalent (i.e. max x is the same in both)
    if isempty(log.scorebins)==1,w.cfmax= 2*10-1; 
    else w.cfmax= 2*size(log.scorebins,1)-1;
    end        
    fcf= @(x,y)(w.cfmax  - abs(x- y) );
    log.define(size(log.define,1)+1, :)={'Conflict' [num2str(w.cfmax) '  - abs diff'] fcf};
    
end


for o=1:1 % Columns for Session 2 
    log.eventratings={'social'; 'vivid';'valence';'arousal';'experience';'happiness';}; 
    col.s2.Trialnum=1;
    col.s2.Option1=2;
    col.s2.Option2=3;
    col.s2.Option=[col.s2.Option1 col.s2.Option2]; 
    col.s2.Option_Onset=4;
    col.s2.ChoiceRT=5;
    col.s2.Motor_Onset=6;  % NOT SURE about this one !!
    col.s2.Choice= 8;   % In excel, 0= Choice 1, 1= Choice 2
    col.s2.PL=[9 10];
    col.s2.PP=[11 12]; 
    col.s2.TrialOK=13; 
    col.s2.Session=14; 
    col.er.social=1;       % Event ratings (row = option)
    col.er.vivid=2; 
    col.er.valence=3; 
    col.er.arousal=4; 
    col.er.experience=5; 
    col.er.happiness=6; 
    
    % Manually calculated by me 
    col.s2_colmax=structmax(col.s2);   
    col.s2.PLPP1= col.s2_colmax+1; 
    col.s2.PLPP2= col.s2_colmax+2; 
    col.s2.PLPP=[col.s2.PLPP1 col.s2.PLPP2];
    col.s2.d1m2.PL=col.s2_colmax+3;  % d1m2: Difference, 1>2 
    col.s2.d1m2.PP=col.s2_colmax+4; 
    col.s2.d1m2.PLPP=col.s2_colmax+5; 
    col.s2.ad1m2.PL=col.s2_colmax+6;  % ad1m2: Absolute Difference, 1>2 
    col.s2.ad1m2.PP=col.s2_colmax+7; 
    col.s2.ad1m2.PLPP=col.s2_colmax+8; 
    col.s2.Trialtype1=col.s2_colmax+9;   % {PL=1 & PP=1}=1, {PL=1 & PP=2}=2
    col.s2.Trialtype2=col.s2_colmax+10; 
    log.Trialtype =reshape(1:100,10,10)';
    col.s2.choPL=col.s2_colmax+11; 
    col.s2.choPP=col.s2_colmax+12; 
    col.s2.choPLPP=col.s2_colmax+13; 
    col.s2.unchoPL=col.s2_colmax+14; 
    col.s2.unchoPP=col.s2_colmax+15; 
    col.s2.unchoPLPP=col.s2_colmax+16;  
    col.s2.conchoPL=col.s2_colmax+17;   % Congruent chioce: 1= PLcho>PLuncho,-1 = opposite, 0=neither
    col.s2.conchoPP=col.s2_colmax+18;   
    col.s2.conchoPLPP=col.s2_colmax+19;   
    col.s2.Trialcf=col.s2_colmax+20;   %  Is this a conflict trial? 
    col.s2.PPbest=col.s2_colmax+21; 
    col.s2.PLbest=col.s2_colmax+22; 
    col.s2.PLPPbest=col.s2_colmax+23; 
    col.s2.PLcf=col.s2_colmax+24;    
    col.s2.PPcf=col.s2_colmax+25;
    col.s2.PLPPcf=col.s2_colmax+26;
    
    
    % Event ratings 
    col.s2_colmax2=structmax(col.s2);   
    col.s2.er.social1 =col.s2_colmax2+1;   
    col.s2.er.social2 =col.s2_colmax2+2;   
    col.s2.er.vivid1 =col.s2_colmax2+3;   
    col.s2.er.vivid2 =col.s2_colmax2+4;   
    col.s2.er.valence1=col.s2_colmax2+5;   
    col.s2.er.valence2=col.s2_colmax2+6;   
    col.s2.er.arousal1=col.s2_colmax2+7;   
    col.s2.er.arousal2=col.s2_colmax2+8;   
    col.s2.er.experience1=col.s2_colmax2+9;   
    col.s2.er.experience2=col.s2_colmax2+10;   
    col.s2.er.happiness1=col.s2_colmax2+11;   
    col.s2.er.happiness2=col.s2_colmax2+12;   
    col.s2.er.d1m2.social=col.s2_colmax2+13;   
    col.s2.er.d1m2.vivid=col.s2_colmax2+14;   
    col.s2.er.d1m2.valence=col.s2_colmax2+15;   
    col.s2.er.d1m2.arousal=col.s2_colmax2+16;   
    col.s2.er.d1m2.experience=col.s2_colmax2+17;   
    col.s2.er.d1m2.happiness=col.s2_colmax2+18;   
    col.s2.er.ad1m2.social=col.s2_colmax2+19;   
    col.s2.er.ad1m2.vivid=col.s2_colmax2+20;   
    col.s2.er.ad1m2.valence=col.s2_colmax2+21;   
    col.s2.er.ad1m2.arousal=col.s2_colmax2+22;   
    col.s2.er.ad1m2.experience=col.s2_colmax2+23;   
    col.s2.er.ad1m2.happiness=col.s2_colmax2+24;
    
end 


if request.do_mainstudy
    
    % [PROGRESS] ################
    %   Note: This script should allow you to recompile the data whenever you
    %   realize you need something else.
    %   Generally, pick ONLY the data that you need from the excel.
    %   Currently only Session 2 (Choice) is compiled (for the modelling). Work
    %           on compiling the others when you get to them.
    %           WHEN YOU DO get to the other sessions, make sure tt
    %            columns specifications are not mixed up between sessions
    %   Columns tt need to be verified w Phoebe are marked in the columns
    %   section!!
    %
    % ############################################
    
    subjdata=[log.subjects  cell(log.n_subjs, 3)];  % Subject, Rating, Choice, ChoiceInLab
    details.readme{1}='Basic data compiled from Phoebe''s excel files (assumed correct). Some diff scores calc''d. SEE SCRIPT do_FormatData.m';
    for s=1: log.n_subjs
        disp(['Subject ' num2str(s) '   -  ' log.subjects{s}  ' ------------------']);
        
        % ##########################################
        % Event ratings  ----------------------------------------
        [ws.n1 ws.t ws.p1]= xlsread([where.beh fs log.subjects{s} fs log.subjects{s} '_behdata.xlsx'],'survey');
        [ws.n2 ws.t ws.p2]= xlsread([where.beh fs log.subjects{s} fs log.subjects{s} '_behdata.xlsx'],'social');
        ws.er= [ws.n2(1:80,3) ws.n1(1:80, 3:7) ];
        
        % ##########################################
        % Session 2 (Choice) ----------------------------------------
        [ws.n ws.t ws.r2]= xlsread([where.beh fs log.subjects{s} fs log.subjects{s} '_behdata.xlsx'],'session2_all');
        ws.rr2=ws.r2(2:end,:); ws.rr2= ws.rr2(:); ws.rr2(find(cellfun(@(x)1- isnumeric(x), ws.rr2)))={nan};   % Nan out blanks
        ws.r2= cell2mat(  reshape(ws.rr2, size(ws.r2)-[1 0])  ); ws.r2(isnan(ws.r2(:,col.s2.Trialnum )), :)=[];
        ws.cols2= [col.s2.Trialnum col.s2.Option1 col.s2.Option2 col.s2.Option_Onset col.s2.ChoiceRT col.s2.Motor_Onset col.s2.Choice col.s2.PL(1)  col.s2.PL(2)  col.s2.PP(1) col.s2.PP(2)];  % Fill in data (only pick what's needed)
        ws.d2(:, ws.cols2) = ws.r2(:, ws.cols2) ;
        ws.d2(:, col.s2.Choice)=ws.d2(:, col.s2.Choice)+1;
        ws.d2(:, col.s2.TrialOK)=  sum([isnan(ws.d2(:, col.s2.Choice)) ws.d2(:, col.s2.Choice)==0  isnan(ws.d2(:, col.s2.PL(1)))    isnan(ws.d2(:, col.s2.PL(2))) isnan(ws.d2(:, col.s2.PP(1)))    isnan(ws.d2(:, col.s2.PP(2)))],2)==0;
        ws.d2(ws.d2(:, col.s2.TrialOK)==1, col.s2.Trialtype1)= cell2mat( cellfun(@(x,y,d)d(x,y),  num2cell(ws.d2(ws.d2(:, col.s2.TrialOK)==1, col.s2.PL(1))) , num2cell(ws.d2(ws.d2(:, col.s2.TrialOK)==1, col.s2.PP(1))), repmat({log.Trialtype},  sum(ws.d2(:, col.s2.TrialOK)==1),1), 'UniformOutput', 0));
        ws.d2(ws.d2(:, col.s2.TrialOK)==1, col.s2.Trialtype2)= cell2mat( cellfun(@(x,y,d)d(x,y),  num2cell(ws.d2(ws.d2(:, col.s2.TrialOK)==1, col.s2.PL(2))) , num2cell(ws.d2(ws.d2(:, col.s2.TrialOK)==1, col.s2.PP(2))), repmat({log.Trialtype},  sum(ws.d2(:, col.s2.TrialOK)==1),1), 'UniformOutput', 0));
        ws.d2(:, col.s2.Session)=2;
        ws.d2(:, col.s2.ChoiceRT) = ws.d2(:, col.s2.Motor_Onset) - ws.r2(:, col.s2.Option_Onset);
        ws.d2(isnan(sum(ws.d2(:,  [col.s2.PL(1) col.s2.PP(1) col.s2.PL(2) col.s2.PP(2)]),2)),  col.s2.TrialOK)=0;
        
        if isempty(log.scorebins)==0
            ws.s2newbins= nan(size(ws.d2,1),4);   % PL1, PL2, PP1, PP2
            for i=1:size(log.scorebins,1)
                ws.binindex= sum(cell2mat( cellfun(@(x, d)d==x, num2cell(log.scorebins{i})', repmat({ws.d2(:, col.s2.PL(1))},  length(log.scorebins{i}),1), 'UniformOutput',0)'),2 );
                ws.s2newbins(find(ws.binindex), 1)=i;
                ws.binindex= sum(cell2mat( cellfun(@(x, d)d==x, num2cell(log.scorebins{i})', repmat({ws.d2(:, col.s2.PL(2))},  length(log.scorebins{i}),1), 'UniformOutput',0)'),2 );
                ws.s2newbins(find(ws.binindex), 2)=i;
                ws.binindex= sum(cell2mat( cellfun(@(x, d)d==x, num2cell(log.scorebins{i})', repmat({ws.d2(:, col.s2.PP(1))},  length(log.scorebins{i}),1), 'UniformOutput',0)'),2 );
                ws.s2newbins(find(ws.binindex), 3)=i;
                ws.binindex= sum(cell2mat( cellfun(@(x, d)d==x, num2cell(log.scorebins{i})', repmat({ws.d2(:, col.s2.PP(2))},  length(log.scorebins{i}),1), 'UniformOutput',0)'),2 );
                ws.s2newbins(find(ws.binindex), 4)=i;
            end
            ws.d2(:, [col.s2.PL(1) col.s2.PL(2) col.s2.PP(1) col.s2.PP(2)])=  ws.s2newbins;
        end
        
        % Calculate new scores from raw
        [ ws.d2] = f_calcfromraw( ws.d2, col.s2);
        ws.d2(:, col.s2.PLcf) = fcf(  ws.d2(:, col.s2.PL(1)), ws.d2(:, col.s2.PL(2)) );
        ws.d2(:, col.s2.PPcf) = fcf(  ws.d2(:, col.s2.PP(1)), ws.d2(:, col.s2.PP(2)) );
        ws.d2(:, col.s2.PLPPcf) = fcf(  ws.d2(:, col.s2.PLPP1), ws.d2(:, col.s2.PLPP2) );

         
        % Event ratings
        ws.d2(ws.d2(:, col.s2.TrialOK)==1, [col.s2.er.social1  col.s2.er.vivid1  col.s2.er.valence1  col.s2.er.arousal1  col.s2.er.experience1  col.s2.er.happiness1]) =  ws.er( ws.d2(ws.d2(:, col.s2.TrialOK)==1, col.s2.Option1), [col.er.social col.er.vivid  col.er.valence col.er.arousal col.er.experience col.er.happiness ]);
        ws.d2(ws.d2(:, col.s2.TrialOK)==1, [col.s2.er.social2 col.s2.er.vivid2 col.s2.er.valence2 col.s2.er.arousal2 col.s2.er.experience2 col.s2.er.happiness2]) = ws.er( ws.d2(ws.d2(:, col.s2.TrialOK)==1, col.s2.Option2), [col.er.social col.er.vivid  col.er.valence col.er.arousal col.er.experience col.er.happiness ]);
        ws.d2(:, [col.s2.er.d1m2.social  col.s2.er.d1m2.vivid  col.s2.er.d1m2.valence  col.s2.er.d1m2.arousal  col.s2.er.d1m2.experience  col.s2.er.d1m2.happiness]) = ws.d2(:,  [col.s2.er.social1  col.s2.er.vivid1  col.s2.er.valence1  col.s2.er.arousal1  col.s2.er.experience1  col.s2.er.happiness1]) - ws.d2(:, [col.s2.er.social2 col.s2.er.vivid2 col.s2.er.valence2 col.s2.er.arousal2 col.s2.er.experience2 col.s2.er.happiness2]);
        ws.d2(:, [col.s2.er.ad1m2.social col.s2.er.ad1m2.vivid col.s2.er.ad1m2.valence col.s2.er.ad1m2.arousal  col.s2.er.ad1m2.experience  col.s2.er.ad1m2.happiness ])= abs(ws.d2(:,  [col.s2.er.d1m2.social  col.s2.er.d1m2.vivid  col.s2.er.d1m2.valence  col.s2.er.d1m2.arousal  col.s2.er.d1m2.experience  col.s2.er.d1m2.happiness]));
        
        
        % ##########################################
        
        
        % Save
        subjdata{s,3}=ws.d2;   s2choice=subjdata{s,3};
        subjdata{s,4}=ws.er;   eventratings=subjdata{s,4};
%         save([where.beh filesep log.subjects{s} filesep log.subjects{s} '_behdata' request.filename_suffix '.mat'], 'col', 's2choice', 'eventratings');
        
        ws=[];
    end
    
    
    
    for o=1:1    % Manually check that values are ok?
        docheck = 0;
        if docheck
            s=10;
            col2=col.s2;
            col2.PL1= col2.PL(1);
            col2.PL2= col2.PL(2);
            col2.PP1= col2.PP(1);
            col2.PP2= col2.PP(2);
            d2= subjdata{s,2+1 };
            col2.ChoiceOrg=   structmax(col2)+1;
            d2(:, col2.ChoiceOrg)= d2(:, col2.Choice);
            d2(:, col2.Choice)= d2(:, col2.Choice)==2; %  -- In excel, 0= Choice 1, 1= Choice 2
            % Session 2, Chosen & unchosen logged properly?
            d2(d2(:, col2.Choice)==0, [col2.choPL col2.choPP col2.unchoPL col2.unchoPP])-d2(d2(:, col2.Choice)==0, [col2.PL1 col2.PP1 col2.PL2 col2.PP2 ])
            
        end
    end
    
    
    % Save overall for modelling
    save([where.where fs '4b Modelling' fs '2 Inputs' fs '2 Data' fs  'All data' request.filename_suffix ' (' date ')'], 'subjdata', 'col', 'log');
    save([where.where fs '4a Beh analysis basic'  fs '2 Data' fs 'All data' request.filename_suffix ' (' date ')'], 'subjdata', 'col', 'log');
    
end



% What are the conflict scores? 
 d{1}=[]; d{2}=[]; 
for s=1:log.n_subjs
    d{1}=[d{1}; unique(subjdata{s,3}(:, [col.s2.PLcf col.s2.PPcf]))]; 
    d{2}=[d{2}; unique(subjdata{s,3}(:, [col.s2.PLPPcf]))]; 
end
d{1}(isnan(d{1}))=[];  d{2}(isnan(d{2}))=[]; 

% unique(d{1}), unique(d{2})





%% Do the same for pilot 

if request.do_pilot
    
    for o=1:1 % Setup
        subjdata=[]; log.datalog= [];         
        log.subjects=[cellfun(@(x)['p0' num2str(x)], num2cell(1:9),'UniformOutput',0)'; cellfun(@(x)['p' num2str(x)], num2cell(10:20),'UniformOutput',0)';];
        log.n_subjs = length(log.subjects);  log.datalog='Pilot data';
        col.s1=  rmfield(col.s2,'er');
        
        % Event ratings
        w=[]; [w.n w.t w.r]= xlsread([where.beh fs   'data_behpilotstudy.xlsx'], 'Rating_ratings');
        d_PL=   w.r(2:81, 3:22)';   d=d_PL; d(cellfun(@(x)ischar(x),  d(: )))= repmat({nan}, sum(cellfun(@(x)ischar(x),  d(: ))),1);
        d_PL=reshape(cell2mat(d), size(d_PL));
        d_PP=  w.r(2:81,27:46)';    d=d_PP; d(cellfun(@(x)ischar(x),  d(: )))= repmat({nan}, sum(cellfun(@(x)ischar(x),  d(: ))),1);
        d_PP=reshape(cell2mat(d), size(d_PP));
        
        % Choices  % d_rt: row=subject, col=trial
        [w.n w.t d_ch]= xlsread([where.beh fs   'data_behpilotstudy.xlsx'], 'Choice_choices');
        [d_rt w.t w.r]= xlsread([where.beh fs   'data_behpilotstudy.xlsx'], 'Choice_rts');  d_rt= d_rt(1:20,:);
        
    end
    
    % s=3;
    % [mean(d_PL(s,:),2) std(d_PL(s,:)')]
    % [mean(d_PP(s,:),2) std(d_PP(s,:)')]
    % PL and PP scores are correct at this point
    
    subjdata= [log.subjects cell(log.n_subjs,1)];
    for s=1:log.n_subjs
        ws.d(:,  col.s1.Trialnum)=1:80;
        ws.ch=  d_ch(3:end, find(strcmp(d_ch(1,:), ['Sub ' num2str(s)])): find(strcmp(d_ch(1,:), ['Sub ' num2str(s)]))+2);  % Opt 1, Opt 2, Choice
        ws.d(:,  [col.s1.Option1 col.s1.Option2 col.s1.Choice])= cell2mat(ws.ch);
        ws.d(:,  col.s1.ChoiceRT)= d_rt(s,2:81)';
        ws.d(:,  col.s1.Option_Onset)=nan;
        ws.d(:,  col.s1.TrialOK)=1;
        ws.d(:,  col.s1.Session)=nan;
        
        for t=1:size(ws.d,1)
            ws.d(t, col.s1.PL(1)) =  d_PL(s, ws.d(t, col.s1.Option1));
            ws.d(t, col.s1.PL(2)) =  d_PL(s, ws.d(t,col.s1.Option2));
            ws.d(t, col.s1.PP(1)) =  d_PP(s, ws.d(t,col.s1.Option1));
            ws.d(t, col.s1.PP(2)) =  d_PP(s, ws.d(t,col.s1.Option2));
        end
        
        %     if sum(isnan(ws.d(:, col.s1.PL(1))))+ sum(isnan(ws.d(:, col.s1.PL(2))))+ sum(isnan(ws.d(:, col.s1.PP(1))))+ sum(isnan(ws.d(:, col.s1.PP(2))))~=0; error; end
            ws.d(isnan(sum(ws.d(:,  [col.s1.PL(1) col.s1.PP(1) col.s1.PL(2) col.s1.PP(2)]),2)),  col.s1.TrialOK)=0;
        
        if isempty(log.scorebins)==0
            ws.s1newbins= nan(size(ws.d,1),4);   % PL1, PL2, PP1, PP2
            for i=1:size(log.scorebins,1)
                ws.binindex= sum(cell2mat( cellfun(@(x, d)d==x, num2cell(log.scorebins{i})', repmat({ws.d(:, col.s1.PL(1))},  length(log.scorebins{i}),1), 'UniformOutput',0)'),2 );
                ws.s1newbins(find(ws.binindex), 1)=i;
                ws.binindex= sum(cell2mat( cellfun(@(x, d)d==x, num2cell(log.scorebins{i})', repmat({ws.d(:, col.s1.PL(2))},  length(log.scorebins{i}),1), 'UniformOutput',0)'),2 );
                ws.s1newbins(find(ws.binindex), 2)=i;
                ws.binindex= sum(cell2mat( cellfun(@(x, d)d==x, num2cell(log.scorebins{i})', repmat({ws.d(:, col.s1.PP(1))},  length(log.scorebins{i}),1), 'UniformOutput',0)'),2 );
                ws.s1newbins(find(ws.binindex), 3)=i;
                ws.binindex= sum(cell2mat( cellfun(@(x, d)d==x, num2cell(log.scorebins{i})', repmat({ws.d(:, col.s1.PP(2))},  length(log.scorebins{i}),1), 'UniformOutput',0)'),2 );
                ws.s1newbins(find(ws.binindex), 4)=i;
            end
            ws.d(:, [col.s1.PL(1) col.s1.PL(2) col.s1.PP(1) col.s1.PP(2)])=  ws.s1newbins;
        end
        
        % Calculate new scores from raw
        [ ws.d] = f_calcfromraw( ws.d, col.s1);
        ws.d(:, col.s1.PLcf) = fcf(  ws.d(:, col.s1.PL(1)), ws.d(:, col.s1.PL(2)) );
        ws.d(:, col.s1.PPcf) = fcf(  ws.d(:, col.s1.PP(1)), ws.d(:, col.s1.PP(2)) );
        ws.d(:, col.s1.PLPPcf) = fcf(  ws.d(:, col.s1.PLPP1), ws.d(:, col.s1.PLPP2) );
        
        subjdata{s,1+1}=ws.d;
        ws=[];
    end
    
    
    % Save overall for modelling
    save([where.where fs '4b Modelling' fs '2 Inputs' fs 'All pilot data' request.filename_suffix ' (' date ')'], 'subjdata', 'col', 'log');
    save([where.where fs '4a Beh analysis basic'  fs '2 Data' fs 'All pilot data' request.filename_suffix ' (' date ')'], 'subjdata', 'col','log');
    
    
end