clc;clear;close all
clear; close all;
User2 = load('/Users/mengxue/Downloads/TmpData/FifteenCm.mat'); 
x = {'U1','U2','U3','U4','U5','U6','U7','U8','U9','U10','U11','U12','U13','U14','U15'};
y = {'U1','U2','U3','U4','U5','U6','U7','U8','U9','U10','U11','U12','U13','U14','U15'};
CM= cell(1,1);
CM{1,1}= User2. A.* 100; %
ConfusionMatrix = CM{1,1};
ConfusionMatrix = roundn(ConfusionMatrix, -2);
h = heatmap(x,y,ConfusionMatrix,'FontSize',20);
 h.FontName = 'Times New Roman';
xlabel('Authenticate result')
ylabel('Ground truth')
set(struct(h).NodeChildren(3), 'XTickLabelRotation', 0);
% save('FifteenConfusionMatrix.mat', CM{1,1}); %Authenticate Result
sum1 = 0;
sumrow = ones(1,15);
for i = 1:15
    sumrow(1,i) = sum( ConfusionMatrix(i,:));
    sum1 = sum1 + ConfusionMatrix(i,i);
end
avg = sum1 / 15;
save('FifteenCm1.mat', 'ConfusionMatrix');