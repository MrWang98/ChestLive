clc;
clear all;
close all;

x1 = [0.0083, 0.0273, 0.0136];
x2 = [0.0083, 0.0273, 0.0136];

figure;
h = boxplot([x1' x2'],'Labels',{'False Accept Rate','False Reject Rate'});
grid on;
set(gca,'FontName','Times New Roman','fontsize',18,'GridLineStyle','--');
set(h,'Linewidth',1.5);