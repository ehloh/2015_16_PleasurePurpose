% Look at res
 clear all; close all hidden, clc
 
 f= load('res_fitmodels_s2 (21-Oct-2015).mat'); 
 g= load('res_gridmodels_s2 (21-Oct-2015)1.mat');
 
%% Load data 
%     g: grid r_res, f: fminunc r_res, 
%     gr/fr: individual subject fits
%     gd/fd: modfit details

modname='bpa02';

gr= g.r_res{strcmp(g.r_res(:,1), modname),2};  fr= f.r_res{strcmp(f.r_res(:,1), modname),2}; 
fd = f.details.models( strcmp(f.details.models(:,1), modname), :); gd = g.details.models( strcmp(f.details.models(:,1), modname), :);
gr(:,1)=1:size(gr,1); fr(:,1)=1:size(fr,1); 

%% Do grid & fit parameters match each other? 

figure('color', 'w')
for p=1:fd{2}
    subplot(1, fd{2}, p)
    scatter(fr(:, 3+p), gr(:, 3+p)), set(gca,'FontSize', 15), lsline
    title(fd{3}{p}, 'FontSize', 20), xlabel('fit'), ylabel('grid') 
end

 
 