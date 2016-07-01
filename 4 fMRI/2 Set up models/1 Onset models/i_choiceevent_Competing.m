function [variables c] = i_choiceevent_Competing(eventtype, data, whichvars, variables, col, c )
% [variables c] = i_rateevent_Competing(data, whichvars, variables, col, c )
% Event type: trial (1), options (2), options til choice (3)
%
% Regressor removal before model estimation looks to remove regressors that
% start with prefix 'o_' 
%
% Session 2/3 events: Fixation (1-3s)  ->   Options (5s)  ->   
%       ->   Choice (until response, max 5s)   ->   Chosen (0.5s)
% ------------------------------------------------------------------------------------------

%% 


eventname={'o_trial', 'o_options', 'o_optionschoice'}; 
eventdur=[col.Duration_Trial col.Duration_Options col.Duration_OptionsChoice]; 
 
for p=1:length(whichvars) 
    variables.names{c}=eventname{eventtype};
    variables.onsets{c}=data(:, col.Onset_Options);
    variables.durations{c}=data(:, eventdur(eventtype));
    
    % Parametric modulators
    variables.pmod(c).name{1}= whichvars{p};
    eval(['variables.pmod(c).param{1}=data(:,col.' whichvars{p} ');'])  % parametric modulators
    variables.pmod(c).poly{1}=1;
    
    %
    c=c+1;
end

end

