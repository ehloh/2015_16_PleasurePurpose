% Make quantities that are needed for fmri analyiss 
clear all;  clc 
cd('/Users/EleanorL/Dropbox/SCRIPPS/7b PLPP fmri/1 Preprocessing')


% i_scanningdetails

nDummyVols=6;
nSlicesPerVol=40; 
TRms=3480/nSlicesPerVol; 

save('i_scanningdetails.mat'); 