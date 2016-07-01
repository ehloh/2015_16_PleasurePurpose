function [ transparval ] = f_transpar(parname, parval, fromto)
%  [ transparval ] = f_transpar(parname, parval, fromto)
%
%     parname   names of parameters (cell)
%     parval   values of parameters (cell; multiple values per param ok)
%     fromto  from/to true parameter vals (string, 'from','to')
%
% --------------------------------------------------------------------------------

% Check inputs
if iscell(parval)==0; parval=num2cell(parval);  outputcell=0; else; outputcell=1; end
if iscell(parname)==0; parname=cellstr(parname); end
if strcmp(fromto, 'from')+strcmp(fromto, 'to')~=1; 
    error('Invalid input fromto (i.e. transform from/to true values). Input ''from'' or ''to'''); 
else; if strcmp(fromto, 'from'); transcol=3; else transcol=2; end
end
if length(parval)~= length(parname); error('No. param vals ~= no. of param names!'); end


% Get parameter defaults
[wm pars wm]=f_modelsettings({'allpars'});
transparval =cell(length(parname),1);

% Transform pars
for p=1:length(parname)
    if sum(strcmp(parname{p}, pars(:,1)))~=1; error(['Could not find requested parameter (' parname{p} '). Check model settings.']); end 
    x=parval{p};
    eval(['transparval{p}=  '  pars{find(strcmp(parname{p}, pars(:,1))), transcol}   ';'])
end

% Reformat to match inputs
transparval=reshape(transparval, size(parval));
if outputcell==0;
    transparval=cell2mat(transparval); 
end



end

