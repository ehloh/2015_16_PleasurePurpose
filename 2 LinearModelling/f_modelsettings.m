function [model_defaults  par_specs models] = f_modelsettings(varargin) 
% [model_defaults  par_details models] = f_modelsettings(whichmodels, parsteps) 
% Generate model settings + fetch requested model details 'models'
% 
% Output variable 'models' lists details of each model -------------------------------------
%           Col 1: Model name
%           Col 2: No. free parameters
%           Col 3: Parameter names
%           Col 4: Constraint transform (constraint implemented in softmax/value fxns)
%           Col 5: Inverse constraint transform transform
%           Col 6: Parameter ranges
%
%   To run an ad-hoc model, add info to model_defaults (name, n free par, par
%   names), if there are new parameters, add into to par_transformations (par,
%   par_transformationsation in, par_transformationsation out)
%   (See function content/comment for details of model_defaults & par_transformations)
%
% --------------------------------------------------------------------------------------------

whichmodels=varargin{1};
if length(varargin)==2; parsteps=varargin{2};
    if parsteps==1; parsteps=2; end
else   parsteps=10; 
end

%% Parameter  specifications
%         Col 1: Parameter name
%         Col 2: Constraint transformations 
%         Col 3: Inverse constraint transformation
%         Col 4: Reasonable parameter range (in true parameter terms)
% 
% Constraint transforms (Col 2) are used to keep parameter values mathematically 
%   kosher, with respect to their usage in the softmax and value functions 
%   (e.g. probabilities between 0 and 1 etc). These transforms are applied within 
%   the softmax/value  functions, and listed here for reference. 
% 
% Inverse constraint transforms (Col 3) are applied if you want to use the existing
%   softmax/value function scripts to generate values/behaviour. Since the scripts
%   themselves apply the transformations, the 'fitted' parameter values must be 
%   inverse transformed first before they enter the softmax/value function scripts
%   themselves (or else, e.g. beta will be exponentiated twice).
% 
% Reasonable parameter range (Col 4), used to (a) seed model fitting iterations, 
%   and (b) specify parameter range in grid search. Because values are in'true'
%   parameter terms, values must be inverse constraint transformed before they
%   are fed (as seed values) into the softmax/value function scripts.
%   Ranges are not binding for the iterative model fits.
  
for o=1:1  % [Archive of of parameter setups]
%     'k'     '0.5+ (1./(2+2.*exp(-x)))'        '- log( 1./(2.*x -1) - 1)'          [0.5+eps 1-eps];      % o: skew softmax to default high pleasure 
    
end

par_specs={       
    'b'     '20./(1+exp(-x))'        '-log( (20./x) - 1)'       [0.01    8];                    % b: softmax beta (beta>0)
    'p'     '1./(2+2.*exp(-x))'      '-log(  1./(2.*x)  -1 )'  [0.001    1/2-0.05];      % p: softmax epsilon ( 0 < epsilon < 1/2; sigmoid function)
    %
	'l'     '1./(2+2.*exp(-x))'      '-log(  1./(2.*x)  -1 )'  [0.001    1/2-0.05];       % l: multiplier of PL scores. Needs to be >0 (weird fits otherwise)
    'a'     '4/(1+exp(-x))-2'        '-log((4./(x+2))-1)'          [-2 2];           % a: weighting of attribute (x some othe variable)
    };


do_vistrans=0;
if do_vistrans  % Visualize transforms 
    for p=1:size(par_specs,1)
        wp.x=linspace(par_specs{p,4}(1),par_specs{p,4}(2),10);
        wp.ffr= inline(par_specs{p,3});
        wp.fto= inline(par_specs{p,2});
        
        % To/From transforms correspond?
        disp([ par_specs{p,1} '   - Transform error: '    num2str( sum(abs(wp.x-wp.fto(wp.ffr(wp.x)))))])
        subplot( ceil(size(par_specs,1)/4), 4, p),
        %     plot(wp.x, wp.ffr(wp.x),'g'); hold on,  plot(wp.ffr(wp.x), wp.fto(wp.ffr(wp.x)),'r')
        plot(wp.x, wp.ffr(wp.x),'g'); hold on,  plot(wp.x, wp.fto(wp.ffr(wp.x)),'r')
        axis square, refline(1,0)
    end
end

 

% Create ranges 
logscalepars={'dummy' }; 
linscalepars = cell2mat(cellfun(@(x)find(1-strcmp(par_specs(:,1), x)), logscalepars, 'UniformOutput',0 )); 
lin_parval = cellfun(@(x)linspace(x(1), x(2), parsteps),  par_specs(linscalepars ,4), 'UniformOutput',0 ); 
logscalepars = cell2mat(cellfun(@(x)find(strcmp(par_specs(:,1), x)), logscalepars, 'UniformOutput',0 )); 
log_parval = cellfun(@(x)logspace( log10(x(1)), log10(x(2)), parsteps),  par_specs(logscalepars,4), 'UniformOutput',0 ); 
par_specs(linscalepars ,4) = lin_parval;
par_specs(logscalepars ,4) = log_parval; 

% The following are NOT free parameters; identity listed here for record
par_fixed={  
    % -------------
    'c'     % conflict  
    'd'     % Absolute difference 
        };
    
    
    
%% Model defaults: Free parameters 
%         Col 1: Model name
%         Col 2: # Free parameters
%         Col 3: Free parametes (in order, according to name)

model_defaults={ 
    'allpars'            []      {'b'; 'p'; 'l'; 'a'};
    
    % -----------------------------------------------
    
    'b01'          [ ]        {'b'};            % Linear sum of PL + PP
    'b02_l'          [ ]        {'b' 'l'}; 
    
    
    'bc01'          [ ]        {'b'};             % Sum cf-weighted PL + PP
    'bc02_l'          [ ]        {'b' 'l'}; 
     
    
    
    
    % -----------------------------------------------
    
    'h01_b_choPL'    []      {'b'}; 
    
    };


% Append lapse (p, pk) models 
model_defaults= [model_defaults;    
    [cellfun(@(x)[x(1) 'p' x(2:end)],  model_defaults(2:end,1), 'UniformOutput',0) cell(size(model_defaults,1)-1,1)  cellfun(@(x)[x(1)  {'p'}  x(2:end)],  model_defaults(2:end,3), 'UniformOutput',0) ]; 
%     [cellfun(@(x)[x(1) 'pk' x(2:end)],  model_defaults(2:end,1), 'UniformOutput',0) cell(size(model_defaults,1)-1,1)  cellfun(@(x)[x(1)  {'p' 'k'}  x(2:end)],  model_defaults(2:end,3), 'UniformOutput',0) ] 
] ;
 

% Add no. pars
model_defaults(:,2)= cellfun( @(x)length(x),  model_defaults(:,3),  'UniformOutput',0);


%% Generate details for requested models: 'models'
%       See function description (above) for description of 'models'

models=cell(length(whichmodels),6); % Read details
for m=1:length(whichmodels)
    if sum(strcmp(model_defaults(:,1), whichmodels{m}))~=0
        wm.modnum=find(strcmp(model_defaults(:,1), whichmodels{m}));
        models{m,1}=model_defaults{wm.modnum,1};
        models{m,2}=model_defaults{wm.modnum,2};
        models{m,3}=model_defaults{wm.modnum,3};
        
        % Fetch details for individual parameters
        for p=1:model_defaults{wm.modnum,2}
            wm.parnum=find(strcmp(par_specs(:,1), model_defaults{wm.modnum,3}{p}));
            models{m,4}{p}=par_specs{wm.parnum,2};
            models{m,5}{p}=par_specs{wm.parnum,3};
            models{m,6}{p}=par_specs{wm.parnum,4};
            
        end
    else error(['Cannot find requested model in f_modelsettings: ' whichmodels{m}])
    end
end


end

