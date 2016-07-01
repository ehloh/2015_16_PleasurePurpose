% Check reg specific subjects 
clear all, clc;


% Request
req.subnum=20; 
req.session=1;   % 1, 2, 3
req.scans=[271 20];  % Probably best to max 3
req.refscans=7;   % 1st scan = 7 (6 dummies)
req.funcprefix = 's8'; 

for o=1:1

where.data ='D:\1 PLPP\1 MRI data\';
req.allsubs ={'p01_YH';'p02_MI';'p03_AY';'p04_AW';'p05_CA';'p06_BT';'p07_HC';'p08_KC';'p09_KJ';'p10_YC';'p11_BI';'p12_AL';'p13_MS';'p14_SK';'p15_MM';'p16_SH';'p17_BB';'p18_WB';'p19_HB';'p20_LZ'};
sub = req.allsubs{req.subnum}; 
    f_mat=@(A,x)A(x); 
end

%%

ws.where = [where.data sub fs '1 Preprocessed' fs 'Func_s' num2str(req.session) filesep]; 

f=spm_select('List', ws.where, ['^' req.funcprefix '.*.00' num2str(req.refscans) '-01.img$']);
matlabbatch{1}.spm.util.checkreg.data{1}= [ws.where f ',1'];
for j=1:length(req.scans)
    f=spm_select('List', ws.where, ['^' req.funcprefix '.*.00' num2str(req.scans(j)) '-01.img$']);
    matlabbatch{1}.spm.util.checkreg.data{j+1}= [ws.where f ',1'];
end
 


spm_jobman('initcfg'); spm_jobman('run' , matlabbatch);