function [nll,pch]=f_nllsoftmax_lapse(x, modelinput, varargin)
% [nll,pch]=f_nllsoftmax_lapse(x, modelinput)
% Calculate neg log-likelihood associated with model
%
%   Inputs
%       x                   Array of parameter points (all unique combinations of parameter values for 
%                               each free parameter) (row= param-point, col=parameter)
%                               [1st param=beta, 2nd param=epsilon]
%
%       modelinput      Value function for model, data, fixed param & column specifications (for data)
%                               i.e. {ModelValueFxn      data   fixedpars     col}
%
%       Z                   Priors for hierarchical fit (Z.mu= mean, Z.nui=variance)
%
%       doprior           Apply penalty for hierarchical fit?
% ------------------------------------------------------------------------------------- 

% Additional non-standard inputs?
if isempty(varargin)==0;
    Z=varargin{1}; doprior=varargin{2};
else doprior=0;
end

% Softmax parameters constrained here - all others constrained in value function
beta=20./(1+exp(-x(:,1)));  % 0 < beta <20
epsilon=1./(2+2.*exp(-x(:,2)));   %  0 < epsilon < 1/2 (sigmoid function)

ValueFxn=modelinput{1}; data=modelinput{2}; fPar=modelinput{3}; col=modelinput{4};

%%

nTrials=size(data,1); nParPoints=size(x,1);

% Apply model's value function to derive values 
%       (V: row=param-point, col=trialnum, 3rd dimension=Option)
xx= x(:, [1:strfind(ValueFxn, 'p')-1  strfind(ValueFxn, 'p')+1:end]);  % Feed it through the corresp valfxn without p
ValueFxn(strfind(ValueFxn, 'p'))=[];  
eval(['V=' ValueFxn '(xx, {[] data, fPar,col});']) 

% Calculate p(Choice) using softmax rule
%           p(Choice)=exp(beta*V(Choice))/(exp(beta*V(Option 1))+exp(beta*V(Option 2))+exp(beta*V(Option 3)));
b2 = repmat(beta, [1 nTrials]);
base=(exp(b2.*squeeze(V(:,:,1)-max(V(:))))  +   exp(b2.*squeeze(V(:,:,2)-max(V(:))))  ); % point in parameter space * trial * option
Vchosen =  repmat( data(:, col.Choice)==1, [1 nParPoints])'.*V(:,:,1) + repmat( data(:, col.Choice)==2, [1 nParPoints])'.*V(:,:,2) ;
pch=exp(b2.*(Vchosen-max(V(:))))./base;

% Episilon: Model accepts some stochasticity in choice
%       p(Choice)=  epsilon + (1-2*epsilon)* p(Choice)      [p(Choice) according to existing softmax calculation]
%       Epsilon is thus bound by chance p(Choice) - in this case, 1/2
pch= repmat(epsilon, [1 nTrials]) + (1-2.*repmat(epsilon, [1 nTrials])).*pch;

% Calculate models negative log-likelihood
nll=  - sum(log(pch),2);


%% [Hierarchical fit] Penalize extreme parameter values 
%   Group parameter distribution (assumed gaussian) is used to constrain individual parameters (by applying a nLL penalty to extreme values)

if doprior
%     keyboard
    % Calculating the penalty, as a function of mean (mu) and variance (nui) of parameter distributions
    LL_penalty=    -1/2 * (x-Z.mu)   *  Z.nui  * (x-Z.mu)' - 1/2*log(2*pi/det(Z.nui));
    nll  = nll  - LL_penalty;
%     
%     LL_penalty=    -1/2 * (x-Z.mu)   *  Z.nui  * (x-Z.mu)' - 1/2*log(2*pi/det(Z.nui));
%     nll  = nll  - LL_penalty;
    
    
    %   Quentin's original: lp = -1/2 * (x-Z.mu)'*Z.nui*(x-Z.mu) - 1/2*log(2*pi/det(Z.nui));
    %                                l  = -l  - sum(lp);
    
%     
%      LL_penalty = log(mvnpdf(x, Z.mu, Z.nui));
%      nll  = nll  - LL_penalty;
end

end


