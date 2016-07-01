function [ LL nPar nStepsPerPar] = f_optimgrid_vect( f, modelinput, ParRanges )
% [ LL nPar nStepsPerPar x] = f_optimgrid_vect( f, modelinput, ParRanges )
%
%   [Input]
%     f                      Name of model (softmax function used to calculate negative log-likelhood)
%     modelinput      Other inputs needed by softmax model (value fxn, & others, fed straight to model)
%                               i.e.: {modelvaluefxn    data   details.fixedpar   col}
%     ParRanges             ParRanges of parameter values to cover, for each free parameter
%                                 Format: Cell vector, full range of for each parameter in each cell
%                                     e.g. {1:10:100     2:20:200     3:30:300} 
%                                     =       Parameter 1 =1,11,21..., Parameter 2=2,22,42..., Parameter 3=3,33,63...
%                            Note: Parameter values will be subject to constraints within the model 
%                                       functions. Values fed in here must thus be inverse transformed first.
%
%   [Output]
%     LL                   log likelihoods associated with each param point for specified model 
%     nPar                   no. of params
%     nStepsPerPar                   
% 
%   Zeb's function to implement grid search (take in vector of parameters)
% --------------------------------------------------------------------------------

% Debug:  m=2; f=@f_nllsoftmax; modelinput={ details.models{m,1} ws.data details.fixedpar col}; ParRanges=details.models{m,6};

nPar = length(ParRanges);
nStepsPerPar=cellfun(@(x)length(x), ParRanges); 
nParPoints = prod(nStepsPerPar); % Total no. of steps (each step = combination of all applicable parameters)

iI=1:nParPoints; % used for nPar = 1:3 only

%% TODO: use a more intelligent chunk size, rather than just chunking the first 3 parameters
%% TODO: test the nPar= 3 and 6 cases

switch nPar
    case 1
        ind=iI';
        x = ParRanges{1}(ind(:,1))';
        [LLt ]= f(x, modelinput);
        LL=LLt ;
    case 2
        [a, b] = ind2sub(nStepsPerPar,iI);
        ind = [a' b'];
        x = [ParRanges{1}(ind(:,1))' ParRanges{2}(ind(:,2))'];
        [LLt ]= f(x, modelinput);
        LL = reshape(LLt, nStepsPerPar(1:2));
    case 3
        [a, b, c] = ind2sub(nStepsPerPar,iI);
        ind = [a' b' c'];
        x = [ParRanges{1}(ind(:,1))' ParRanges{2}(ind(:,2))' ParRanges{3}(ind(:,3))'];
        [LLt ]= f(x, modelinput);
        LL = reshape(LLt, nStepsPerPar(1:3));
    case 4
        for iP4=1:nStepsPerPar(4)
            Lt = prod(nStepsPerPar(1:3));
            [a, b, c] = ind2sub(nStepsPerPar(1:3),1:Lt);
            ind = [a' b' c'];
            x = [ParRanges{1}(ind(:,1))' ParRanges{2}(ind(:,2))' ParRanges{3}(ind(:,3))' repmat(ParRanges{4}(iP4), [Lt 1])];
            LLt = f(x, modelinput);
            LL(:,:,:,iP4) = reshape(LLt, nStepsPerPar(1:3));
        end
    case 5
        for iP4=1:nStepsPerPar(4)
            for iP5=1:nStepsPerPar(5)
                Lt = prod(nStepsPerPar(1:3));
                [a, b, c] = ind2sub(nStepsPerPar(1:3),1:Lt);
                ind = [a' b' c'];
                x = [ParRanges{1}(ind(:,1))' ParRanges{2}(ind(:,2))' ParRanges{3}(ind(:,3))' repmat(ParRanges{4}(iP4), [Lt 1]) repmat(ParRanges{5}(iP5), [Lt 1])];
                LLt = f(x, modelinput);
                LL(:,:,:,iP4,iP5) = reshape(LLt, nStepsPerPar(1:3));
            end
        end
    case 6
        for iP4=1:nStepsPerPar(4)
            for iP5=1:nStepsPerPar(5)
                for iP6=1:nStepsPerPar(6)
                    Lt = prod(nStepsPerPar(1:3));
                    [a, b, c] = ind2sub(nStepsPerPar(1:3),1:Lt);
                    ind = [a' b' c'];
                    x = [ParRanges{1}(ind(:,1))' ParRanges{2}(ind(:,2))' ParRanges{3}(ind(:,3))' repmat(ParRanges{4}(iP4), [Lt 1]) repmat(ParRanges{5}(iP5), [Lt 1]) repmat(ParRanges{6}(iP6), [Lt 1])];
                    LLt = f(x, modelinput);
                    LL(:,:,:,iP4,iP5,iP6) = reshape(LLt, nStepsPerPar(1:3));
                end
            end
        end
    case 7        
        for iP4=1:nStepsPerPar(4)
            for iP5=1:nStepsPerPar(5)
                for iP6=1:nStepsPerPar(6)
                    for iP7=1:nStepsPerPar(7)
                        Lt = prod(nStepsPerPar(1:3));
                        [a, b, c] = ind2sub(nStepsPerPar(1:3),1:Lt);
                        ind = [a' b' c'];
                        %
                        x = [ParRanges{1}(ind(:,1))' ParRanges{2}(ind(:,2))' ParRanges{3}(ind(:,3))' repmat(ParRanges{4}(iP4), [Lt 1]) repmat(ParRanges{5}(iP5), [Lt 1]) repmat(ParRanges{6}(iP6), [Lt 1])   repmat(ParRanges{7}(iP7), [Lt 1])   ];
                        LLt = f(x, modelinput);
                        LL(:,:,:,iP4,iP5,iP6,iP7) = reshape(LLt, nStepsPerPar(1:3));
                    end
                end
            end
        end
    otherwise
        error('Invalid no. of params (1-7 only)')
end

end