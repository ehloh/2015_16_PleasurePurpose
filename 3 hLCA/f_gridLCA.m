function [  LL nPar nStepsPerPar] = f_gridLCA( f, modelinput, ParRanges )
%   [Input]
%     f                      Name of model 
%     modelinput      Other inputs needed by softmax model (fed straight to model) 
%     ParRanges             ParRanges of parameter values to cover, for each free parameter
%                                 Format: Cell vector, full range of for each parameter in each cell
%                                     e.g. {1:10:100     2:20:200     3:30:300} 
%                                     =       Parameter 1 =1,11,21..., Parameter 2=2,22,42..., Parameter 3=3,33,63...
%                            Note: Parameter values will be subject to constraints within the model 
%                                       functions. Values fed in here must thus be inverse transformed first.
%   [Output]
%     LL                   log likelihoods associated with each param point for specified model 
%     nPar                   no. of params
%     nStepsPerPar                   
% 
% ---------------------------------------------------------------------------------------------------------
 
nPar= length(ParRanges); 
nStepsPerPar=cellfun(@(x)length(x), ParRanges); 
if power(nStepsPerPar(1), nPar)~=prod(nStepsPerPar); error('Grid script assumes equal parsteps for each parameter, but script is not set up for this'); end 
fmat= @(mat, index)mat(index);
eval(['ParComb= allcomb('  fmat( repmat([num2str(1) ':' num2str(nStepsPerPar(1)) ','], 1, nPar),  1:length(repmat([num2str(1) ':' num2str(nStepsPerPar(1)) ','], 1, nPar))-1) ');'])

% Run model through every point in parameter space 
eval(['LL=nan(' [strrep(strjoin(cellfun(@(x)[num2str(x) ','], num2cell(nStepsPerPar(1:end-1)), 'UniformOutput',0)),',,',',') num2str(nStepsPerPar(end))] ');'])   % This is so long because strjoin works differently on PC vs Mac -__-
for iPoint=1:size(ParComb,1)
    pointL=  f(cellfun(@(ra,i)ra(i), ParRanges, num2cell(ParComb(iPoint, :))), modelinput);   
    eval(['LL('  [strrep(strjoin(cellfun(@(x)[num2str(x) ','], num2cell(ParComb(iPoint,1:size(ParComb,2)-1)), 'UniformOutput',0)),',,',',') num2str(ParComb(iPoint,end))]         ')=pointL;'])
     
    % Display progress 
    if sum(iPoint == 0:50:size(ParComb,1)),  disp(['     grid progress: ' num2str(iPoint) ' out of '  num2str(size(ParComb,1)) ' points in param space  (' num2str(fmat(clock,4)) ':' num2str(fmat(clock,5)) ' hrs)']); end    
end 
end

