function [variables c] = i_mainevents(sessiontype, data, variables, col,c ) 
% Main events (Rating or Choice session) 
 
switch sessiontype
    % Session 1 events --------------------------------------------------
    case 1
        variables.names{c}='event';
        variables.onsets{c}=data(:, col.Onset_event);
        variables.durations{c}= data(:, col.Duration_event);
        
    % Session 2/3 events --------------------------------------------------
    case 2.1  % Entire trial epoch (5s considering options, 0-5s making choice, 0.5s see chosen)
        variables.names{c}='trial';
        variables.onsets{c}=data(:, col.Onset_Options);
        variables.durations{c}=data(:, col.Duration_Trial); 
    case 2.2  % Considering options only 
        variables.names{c}='options';
        variables.onsets{c}=data(:, col.Onset_Options);
        variables.durations{c}=data(:, col.Duration_Options);
    case 2.3 % Considering options + making choice  
        variables.names{c}='optionschoice';
        variables.onsets{c}=data(:, col.Onset_Options);
        variables.durations{c}=data(:, col.Duration_OptionsChoice);  
    otherwise, error('Type of main event not set up yet!');
end 
c=c+1;  

end

 