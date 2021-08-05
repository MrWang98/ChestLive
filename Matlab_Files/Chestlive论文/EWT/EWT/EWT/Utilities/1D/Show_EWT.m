function Show_EWT(ewt,f,rec,EWT_Tic)

% ====================================================================
% function Show_EWT(ewt,f,rec)
%
% This function plots the successive filtered components (low scale 
% first and then wavelets scales). The original and
% reconstructed signals are plotted on a different graph.
% If f and rec are provided, it also plot the original and reconstructed
% signals on a separate figure
%
% Inputs:
%   -ewt: EWT components
%   -f: input signal  (OPTIONNAL)
%   -rec: reconstructed signal (OPTIONNAL)
%
% Author: Jerome Gilles
% Institution: UCLA - Department of Mathematics
% Year: 2012
% Version: 1.0
% =====================================================================

%These lines plot the EWT components
figure;
% x=0:1/length(ewt{1}):(length(ewt{1})-1)/length(ewt{1});
x = EWT_Tic;
l=1;
if length(ewt)>6 %raw is 6
    lm=6; % raw is 6
else
    lm=length(ewt);
end

for k=1:length(ewt)
%    hold on; subplot(lm,1,l); plot(x,ewt{k}); %axis off;
%    if mod(k,6) == 0 %raw is 6
%         figure;
%         l=1;
%    else
%     l=l+1;
%    end
    figure(1234567);
    subplot(4,1,1);
    plot(x,f, 'LineWidth',2);
    ylabel("Origin",'fontsize',22,'FontName','Times New Roman') ;
    set(gca,'linewidth',1.5,'FontName','Times New Roman','fontsize',18,'GridLineStyle','--');
    subplot(4,1,2);
    plot(x,ewt{1}, 'LineWidth',2);
    ylabel('Ewt-C1','fontsize',22,'FontName','Times New Roman') ;
    set(gca,'linewidth',1.5,'FontName','Times New Roman','fontsize',18,'GridLineStyle','--');
    subplot(4,1,3);
    plot(x,ewt{2}, 'LineWidth',2);
    ylabel('Ewt-C2','fontsize',22,'FontName','Times New Roman') ;
% set(gca,'YTick', [-0.005:0.02:0.005]);
%     axes('YTickLabel',{'-0.005','0.005'});
    set(gca,'linewidth',1.5,'FontName','Times New Roman','fontsize',18,'GridLineStyle','--');
    subplot(4,1,4);
    plot(x,ewt{3}, 'LineWidth',2);
    ylabel('Ewt-C3','fontsize',22,'FontName','Times New Roman') ;
    xlabel("Time (s)",'fontsize',22,'FontName','Times New Roman') ;
    set(gca,'linewidth',1.5,'FontName','Times New Roman','fontsize',18,'GridLineStyle','--');
end

%These lines plot f and its reconstruction
if nargin>1
    figure;
    subplot(2,1,1);plot(x,f);
    subplot(2,1,2);plot(x,rec);
end
