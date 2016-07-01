function [variables c] = modname(datastruc, colstruc, c) 
% Non-orthogonalized: Each variable has its own onset + pmod, Duplicate onsets are ignored in first-level contrasts
% 
% datastruc: {1} Rating {2} Choice {3} Choice Lab {4} Choice All
%
% Session 1 events: Fixation (0.5-3s)  ->   Event (5s)  ->   
%       ->   Pleasure rating (until response, max 5s)  ->   Jitter (0.5-3s)  ->   
%       ->   Purpose rating (until response, max 5s)
%  -----------------------------------------------------------------------------

%% 

variables = []; 

% Events presented for rating + pmods (competing)
[variables c] = i_mainevents(1, datastruc{1}(datastruc{1}(:, colstruc{1}.TrialOK)==1,:),  variables, colstruc{1}, c );  

% Events presented for rating + pmods (competing)
[variables c] = i_rateevent_Competing(datastruc{1}(datastruc{1}(:, colstruc{1}.TrialOK)==1,:), {'PLbinrating','PPbinrating'}, variables, colstruc{1}, c );  

% NO-INTEREST regressors
[variables c] = i_nointerest(1, datastruc{1}, {'PLrate'; 'PPrate'; 'Motor'; 'PLrateMiss'; 'PPrateMiss'}, variables, colstruc{1}, c);

   
end

