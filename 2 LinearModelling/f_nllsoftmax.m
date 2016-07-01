function [nll,pch]=f_nllsoftmax(x, modelinput, varargin)
% [nll,pch]=f_nllsoftmax(x, modelinput, varargin)
% Calculate neg log-likelihood associated with model
%
%   Inputs
%       x                   Array of parameter points (all unique combinations of parameter values for 
%                               each free parameter) (row= param-point, col=parameter)
%
%       modelinput    Value function for model, data, fixed param & column specifications (for data)
%                               i.e. {ModelValueFxn      data   fixedpars     col}
%
%       Z                   Priors for hierarchical fit (Z.mu= mean, Z.nui=variance)
%
%       doprior           Apply penalty for hierarchical fit?
% ------------------------------------------------------------------------------------- 

% Additional non-standard inputs?
if isempty(varargin)==0; Z=varargin{1};doprior=varargin{2};
else doprior=0;
end

% beta=exp(x(:,1)); % Beta is constrained here - all other parameters constrained in value function
beta=20./(1+exp(-x(:,1)));
ValueFxn=modelinput{1}; data=modelinput{2}; fPar=modelinput{3}; col=modelinput{4};

%%

nTrials=size(data,1); nParPoints=size(x,1);

% Apply model's value function to derive values 
%       (V: row=param-point, col=trialnum, 3rd dimension=Option)
eval(['V=' ValueFxn '(x, {[] data, fPar,col});']) 

% Calculate p(Choice) using softmax rule
%           p(Choice)=exp(beta*V(Choice))/(exp(beta*V(Option 1))+exp(beta*V(Option 2))+exp(beta*V(Option 3)));
%           Values are shifted by max to avoid inf exp(V)
b2 = repmat(beta, [1 nTrials]);
base=(exp(b2.*squeeze(V(:,:,1)-max(V(:))))  +   exp(b2.*squeeze(V(:,:,2)-max(V(:))))  ); % point in parameter space * trial * option
Vchosen =  repmat( data(:, col.Choice)==1, [1 nParPoints])'.*V(:,:,1) + repmat( data(:, col.Choice)==2, [1 nParPoints])'.*V(:,:,2) ;
pch=exp(b2.*(Vchosen-max(V(:))))./base;

% Calculate models negative log-likelihood
nll=  - sum(log(pch),2);

%% [Hierarchical fit] Penalize extreme parameter values 
%   Group parameter distribution (assumed gaussian) is used to constrain individual parameters (by applying a nLL penalty to extreme values)
%   Make sure Z.mu is in the same space as x

if doprior
    % Calculating the penalty, as a function of mean (mu) and variance (nui) of parameter distributions
    LL_penalty=    -1/2 * (x-Z.mu)   *  Z.nui  * (x-Z.mu)' - 1/2*log(2*pi/det(Z.nui));
    nll = nll  - LL_penalty;
    % CHECK THAT THIS FORMULATION IS CORRECT !!
end

end


