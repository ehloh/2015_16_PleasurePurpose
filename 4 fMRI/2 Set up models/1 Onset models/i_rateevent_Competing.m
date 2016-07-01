function [variables c] = i_rateevent_Competing(data, whichvars, variables, col, c )
% [variables c] = i_rateevent_Competing(data, whichvars, variables, col, c )
% Prefix for onset events: oe_
%
% Session 1 events: Fixation (0.5-3s)  ->   Event (5s)  ->   
%       ->   Pleasure rating (until response, max 5s)  ->   Jitter (0.5-3s)  ->   
%       ->   Purpose rating (until response, max 5s)
% ------------------------------------------------------------------------------------------

%% 

for p=1:length(whichvars)
    variables.names{c}=['o_event'];
    variables.onsets{c}=data(:,col.Onset_event); 
    variables.durations{c}=data(:,col.Duration_event); 
    
    % Parametric modulators
    variables.pmod(c).name{1}= whichvars{p};
    eval(['variables.pmod(c).param{1}=data(:,col.' whichvars{p} ');'])  % parametric modulators
    variables.pmod(c).poly{1}=1;
    
    %
    c=c+1;
end

end

