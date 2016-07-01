% Look at likelihood surface (subject grid) 
% clear all;close all hidden; clc;   
clear all; clc;   
  
% Request 
req.file='p10_YC_se2stp8sam300_modlta01_fit1.mat'; 

for o1=1:1 % Set up  
    w=pwd; if strcmp(w(1),'/'), where.where='/Users/EleanorL/Dropbox/SCRIPPS/3 Pleasure purpose';    else  where.where='C:\Users\e.loh\Dropbox\SCRIPPS\3 Pleasure purpose'; end 
    where.mod =[where.where filesep '4b Modelling']; where.beh=[where.where filesep '3 Behaviour']; 
    path(pathdef), addpath(where.where), cd(where.mod)
    addpath([where.mod filesep '1 Value functions'])  
    addpath([where.mod filesep '1 Value functions' filesep 'hLCA'])  
    
    %
     fmat= @(mat, index)mat(index);  
     
    % Load + parse details 
    req.f= load([where.mod filesep '2 Inputs' filesep '3 Subgrids' filesep req.file]);
    logg.model=  req.file(strfind(req.file, '_mod')+4: strfind(req.file, '_fit')-1); 
    logg.modset= req.f.details.models( strcmp(req.f.details.models(:,1), logg.model),:);
    logg.fp= req.f.details.fixedpar;
    logg.par = [logg.modset{3}'  logg.modset{6}'  logg.modset{7}'];  % parname, parstepvals (true param values), parstepvals (reverse transformed)
    logg.par_trans =  [logg.modset{4}' logg.modset{5}'];
    logg.n_parsteps =  req.f.details.n_parsteps;  logg.n_par = length(logg.par(:,1));
    logg.n_samples = req.f.details.fixedpar.nsamples; 
    logg.n_trials = req.f.details.subj_ntrials;
    nl=  req.f.L;  
    
    % Parse file details 
    logg.sub=  req.file(1:6); 
    
end


%% View  grids for each subject
 
% ind(ip1,ip2.. ipn)  = indices of best point in param space 
try ind = subfit.minL_indices; catch;  eval([  '[' strrep(strjoin(cellfun(@(x)[' ind(' num2str(x) ')'],num2cell(1:logg.n_par), 'UniformOutput',0)),',,',',')  ']= ind2sub(size(nl), find(nl==min(nl(:))));']); end      
  
 
%  Plot likelihood surface 
chance =  -(log(0.5)*logg.n_trials ); 
figure, surf(squeeze(nl(:,ind(2), :)))
title([logg.sub  '     '   logg.model '   (' num2str(logg.n_samples) ' samples)'], 'FontSize',15), xlabel('Dim B', 'FontSize',15), ylabel('Dim A', 'FontSize',15), zlabel( ['L  (chance: ' num2str(chance) ')'], 'FontSize',15) 

 

