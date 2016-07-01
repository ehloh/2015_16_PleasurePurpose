function [variables c] = i_nointerest (sessiontype, data, NoInterestVar, variables, col,c ) 
% No interest regressors for rating/choice sessions (specified)

%% Format data  

switch sessiontype
    case 1 
        d_PLok= data(~isnan(data(:, col.PL)),:);
        d_PPok= data(~isnan(data(:, col.PP)),:);
        d.PLrate= d_PLok(:, col.Onset_PLrate);    % VAS rating events
        d.PPrate= d_PPok(:, col.Onset_PPrate);
        d.Motor=  sortrows([data(:,  col.Onset_PLmotor); data(:,  col.Onset_PPmotor)]);
        d.PLrateMiss = data(isnan(data(:, col.PL)) , col.Onset_PLrate);         % Errors/Omissions
        d.PPrateMiss = data(isnan(data(:, col.PP)) , col.Onset_PPrate);
        if isempty(d.PLrateMiss) && sum(strcmp( NoInterestVar, 'PLrateMiss')) >0, NoInterestVar= NoInterestVar(~strcmp( NoInterestVar, 'PLrateMiss')); end 
        if isempty(d.PPrateMiss) && sum(strcmp( NoInterestVar, 'PPrateMiss')) >0, NoInterestVar= NoInterestVar(~strcmp( NoInterestVar, 'PPrateMiss')); end 
        RegsWSpecialDuration ={}; 
    case 2  
        d.Motor = data(~isnan(data(:, col.Onset_Choice)), col.Onset_Choice);  
        d.ChoiceMiss=  data(isnan(data(:, col.Choice)), col.Onset_Options); 
        d.RatingMissingTrial = data( find(ceil( isnan(data(:, col.PL1))+ isnan(data(:, col.PL2))+ isnan(data(:, col.PP1)) + isnan(data(:, col.PP2)))), col.Onset_Options);
        d.RatingMissingTrial_duration = data( find(ceil( isnan(data(:, col.PL1))+ isnan(data(:, col.PL2))+ isnan(data(:, col.PP1)) + isnan(data(:, col.PP2)))), col.Duration_Trial);
        
        if isempty(d.RatingMissingTrial) && sum(strcmp( NoInterestVar, 'RatingMissingTrial')) >0, NoInterestVar= NoInterestVar(~strcmp( NoInterestVar, 'RatingMissingTrial')); end
        if isempty(d.ChoiceMiss) && sum(strcmp( NoInterestVar, 'ChoiceMiss')) >0, NoInterestVar= NoInterestVar(~strcmp( NoInterestVar, 'ChoiceMiss')); end 
        
        
        % Regressors where duration must be read out manually (*_duration)
        RegsWSpecialDuration ={'RatingMissTrial'; 'test'}; 
    otherwise, error('Session not defined!')
        
end 

%% Construct no interest regressors   

for p=1:length(NoInterestVar) 
    variables.names{c}=['n_' NoInterestVar{p}];
    eval(['variables.onsets{c}=d.' NoInterestVar{p} ';'])
     
    % Non-zero duration?
    if sum(strcmp( RegsWSpecialDuration, NoInterestVar{p}))>0,  eval(['variables.durations{c}=d.' NoInterestVar{p} '_duration;'])
    else variables.durations{c}=0;
    end
    
    % Parametric modulators (for particular no-interest variables)
    if strcmp(NoInterestVar{p}, 'OutcomePresented')==1
        disp('Not set up for pmods yet!')
        variables.pmod(c).name{1}='n_OutcomePresentedMagnitude';
        variables.pmod(c).param{1}=d.OutcomePresentedMagnitude;
        variables.pmod(c).poly{1}=1;
    end
    %
    c=c+1;
end


end

