% Fetch all data for 1 subject to check if correct
clear all



ons= load('/Users/EleanorL/Dropbox/WorkPC/PLPP/p01_YH/2 First level/p01_YH_onsets_r1_PLnPPrating.mat');
spm= load('/Users/EleanorL/Dropbox/WorkPC/PLPP/p01_YH/2 First level/r1_PLnPPrating Contrasted/SPM.mat');spm=spm.SPM;
b= load('/Users/EleanorL/Dropbox/SCRIPPS/3 Pleasure purpose/3 Behaviour/p01_YH/p01_YH_behdata.mat');


% Onsets? 

n=spm.xX.name'; openvar n

spm

