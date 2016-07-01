function [ out ] = f_subsample(varargin)
% [ out ] = f_subsample(varargin)
%  Full inputs (2): Fetch list of subjects (cell) 
%         sublist= f_subsample(FLmodel, SortCriteria);
%                       e.g.  f_subsample('r1b', 'oe_PLvPP', 1) 
%                       Order is always from low/negative to high r values 
%
%  No inputs (0): Output list (cell) of models and their available SortCriteria 
%         {FLmod SortCriteria;...} =  f_subsample( ) 
%                                               e.g.  {'r1' 'oe_PLvPP'; 'r1' 'PLvPP'; ..}
%
%  1 input (FLmodel):  Output available SortCriteria for requested model  
%         {SortCriteria} =  f_subsample(FLmod) 
%                                e.g. {'oe_PLvPP';'PLvPP'  ..} = f_subsample('r1b') 
%
%
%  Note that 'beh' is also an option (in place of FLmodel)
%  FLmod can also be input as entire FL name 
% --------------------------------------------------------------------------------

switch length(varargin)
    case 0, do='ListFLxSortCrit';
    case 1, do='ListSortCrit'; mod = varargin{1}; 
    case 2, do='ListSub'; 
        mod = varargin{1}; 
        sortcrit= varargin{2};  
    otherwise; error('Invalid n inputs for fetching subsamples') 
end 
f_mat=@(A,x)A(x); 
          
%% Subjects subsamples: Always listed from low/negative correlation to high/positive correlation 
%   It only makes sense to apply subsamples for criteria where subjects show a lot of spread 
%   See script s_CheckRegCorr to generate these subsamples 


% Behaviour session 1 (Rating)
subs.beh1.PLvPP= {'p04_AW';'p09_KJ';'p12_AL';'p01_YH';'p08_KC';'p11_BI';'p17_BB';'p18_WB';'p19_HB';'p02_MI';'p14_SK';'p15_MM';'p06_BT';'p10_YC';'p13_MS';'p16_SH';'p03_AY';'p07_HC';'p05_CA';'p20_LZ';};
subs.beh1.PLvPP_abs= {'p19_HB';'p18_WB';'p02_MI';'p17_BB';'p11_BI';'p08_KC';'p14_SK';'p01_YH';'p12_AL';'p09_KJ';'p04_AW';'p15_MM';'p06_BT';'p10_YC';'p13_MS';'p16_SH';'p03_AY';'p07_HC';'p05_CA';'p20_LZ';}; 
subs.beh1.PLvPP_abs_r04= {'p19_HB';'p18_WB';'p02_MI';'p17_BB';'p11_BI';'p08_KC';'p14_SK';'p01_YH';'p12_AL';'p09_KJ';'p04_AW';'p15_MM';}; 

% Behaviour session 2 (Choice) 


% r1 First Level 
subs.r1.regr_oe_PLvPP=  {'p04_AW';'p12_AL';'p09_KJ';'p01_YH';'p08_KC';'p11_BI';'p17_BB';'p18_WB';'p19_HB';'p02_MI';'p14_SK';'p15_MM';'p06_BT';'p10_YC';'p13_MS';'p16_SH';'p03_AY';'p07_HC';'p05_CA';'p20_LZ';}; 
subs.r1.regr_oe_PLvPP_abs= {'p19_HB';'p18_WB';'p02_MI';'p17_BB';'p11_BI';'p08_KC';'p14_SK';'p01_YH';'p09_KJ';'p12_AL';'p04_AW';'p15_MM';'p06_BT';'p10_YC';'p13_MS';'p16_SH';'p03_AY';'p07_HC';'p05_CA';'p20_LZ';};
subs.r1.regr_oe_PLvPP_abs_r04= {'p19_HB';'p18_WB';'p02_MI';'p17_BB';'p11_BI';'p08_KC';'p14_SK';'p01_YH';'p09_KJ';'p12_AL';'p04_AW';'p15_MM';};  % absolute correlation r <0.4. n=12 

% c2 First level  
subs.c2.regr_ot_PLvPP_abs ={'p18_WB';'p17_BB';'p01_YH';'p14_SK';'p08_KC';'p04_AW';'p02_MI';'p19_HB';'p11_BI';'p09_KJ';'p12_AL';'p06_BT';'p15_MM';'p03_AY';'p05_CA';'p16_SH';'p07_HC';'p10_YC';'p13_MS';'p20_LZ';}; 
subs.c2.regr_ot_PLvPP_abs_r04 ={'p18_WB';'p17_BB';'p01_YH';'p14_SK';'p08_KC';'p04_AW';'p02_MI';'p19_HB';'p11_BI';'p09_KJ';'p12_AL';'p06_BT';}; 
subs.c2b.regr_ot_PLvPP_abs= {'p11_BI';'p18_WB';'p08_KC';'p17_BB';'p14_SK';'p12_AL';'p02_MI';'p04_AW';'p01_YH';'p19_HB';'p09_KJ';'p05_CA';'p06_BT';'p03_AY';'p15_MM';'p10_YC';'p16_SH';'p13_MS';'p07_HC';'p20_LZ';};
subs.c2b.regr_ot_PLvPP_abs_r04= {'p11_BI';'p18_WB';'p08_KC';'p17_BB';'p14_SK';'p12_AL';'p02_MI';'p04_AW';'p01_YH';'p19_HB';'p09_KJ';'p05_CA';'p06_BT';}; 


 
% subs.c3.regr_ot_PLcfvPPcf_abs = 
%% Task: ListFLxSortCrit
% List all FLs and their associated sorting criteria 

if strcmp(do, 'ListFLxSortCrit')
    modcrit={};  mods = fieldnames(subs);
    for m=1:length(mods) 
        crits =  getfield(subs, mods{m});  
        modcrit=[modcrit; [repmat(mods(m), length(fieldnames(crits)),1) fieldnames(crits)]] ; 
    end 
    out=modcrit;
end 

%% Task: SortCrit
% List sorting criteria for requested FL 

if strcmp(do,  'ListSortCrit')
    if ~isempty(strfind(mod, '_')), mod=  mod(1:f_mat(strfind(mod, '_'),1)-1); end  % Shorten name if full FL model name was input
    if isfield(subs, mod)==0; error('Could not find requested model'); end
    crits = fieldnames(subs.(mod));
    out= crits;
end 

%% Task: ListSub
% List subjects for requested FL model & sort criteria 

if strcmp(do, 'ListSub')
    if ~isempty(strfind(mod, '_')), mod=  mod(1:f_mat(strfind(mod, '_'),1)-1); end  % Shorten name if full FL model name was input
    if strcmp(sortcrit(1:3), 'beh'), mod=['beh' sortcrit(4)];  sortcrit = sortcrit(strfind(sortcrit, '.')+1:end); end
    
    if isfield(subs, mod)==0; error('Could not find requested model'); end
    if sum(strcmp(fieldnames(subs.(mod)),sortcrit))==0; error('Could not find requested sorting criteria (for requested FL model)'); end
    sublist= subs.(mod).(sortcrit); 
    out= sublist;  
end 
 

end

