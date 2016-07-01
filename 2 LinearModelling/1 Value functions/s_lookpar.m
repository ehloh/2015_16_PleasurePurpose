%% Work out parameters 
clear all, clc

%%  i parameter 

p=0.2; 
scavals = [0 0.25 0.5 0.75 1];  % i parameter value 
x=-10:0.05:10;
sf=@(x)1./(1+exp(-x)); 
sfp = @(x,p)p + (1-2*p).*sf(x);
sfp_scale = @(x,p, sc) 2*p*sc + (1-2*p).*sf(x);

close all, figure('color','w')
% subplot(2,2,1:2)  % General formulation 
plot(x, ones(length(x),1)*0.5, 'color',[0 0 0]); hold on, plot([0 0 0], [0 0.5 1], 'color',[0 0 0]); hold on
for i=1:length(scavals)
    plot(x, sfp_scale(x,p, scavals(i)), 'color', [ i/5 0 0], 'linewidth', 2); hold on
    text(5.6, 0.1+i*0.05, ['K=' num2str(scavals(i))], 'color',[ i/5 0 0],'fontsize',20)
end
text(-9.5, 0.95,  ['p=' num2str(p)] , 'color', [0 0 0],'fontsize',20); title('K parameter','fontsize',30)
set(gca,'fontsize',20),  xlabel('v(A)-v(B)'), ylabel('p(A)')

% 
% 
% subplot(2,2,3)  % General formulation 
% plot(x, ones(length(x),1)*0.5, 'color',[0 0 0]); hold on; plot([0 0 0], [0 0.5 1], 'color',[0 0 0]); hold on
% for i=1:length(scavals)
%     plot(x, sfp_scale(x,p, scavals(i)), 'color', [ i/5 0 0], 'linewidth', 2); hold on
%     text(5.6, 0.1+i*0.05, ['i=' num2str(scavals(i))], 'color',[ i/5 0 0],'fontsize',20)
% end
% 
% text(-9.5, 0.95,  ['p=' num2str(p)] , 'color', [0 0 0],'fontsize',20)
% set(gca,'fontsize',20),  xlabel('v(PL)-v(B)'), ylabel('p(A)')
% title('PL','fontsize',30)






