% Explain the data positon. and the data in the file.
%bins[1850]	bins[1851]	bins[1852]	bins[1853]
%bins[1854]	bins[1855]	bins[1856]	bins[1857]
%bins[1858]	bins[1859]	bins[1860]	bins[1861]
%bins[1862]	bins[1863]	bins[1864]	bins[1865]
% figure;plot(xtic,bin1650);figure;plot(xtic,bin1651);figure;plot(xtic,bin1652);figure;plot(xtic,bin1653);
% figure;plot(xtic,bin1654);figure;plot(xtic,bin1655);figure;plot(xtic,bin1656);figure;plot(xtic,bin1657);
% figure;plot(xtic,bin1658);figure;plot(xtic,bin1659);figure;plot(xtic,bin1660);figure;plot(xtic,bin1661);
clear,clc;
close all;

% Read data from files (1)/Users/mengxue/Downloads/
vM = csvread('/Users/mengxue/Downloads/TEST_BINS(2).csv', 2,0);
time = vM(:,1);
xtic= (time(:) - time(1))./1000;
bin1650 = vM(:,2);  bin1651 = vM(:,3);  bin1652 = vM(:,4);  bin1653 = vM(:,5);
bin1654 = vM(:,6);  bin1655 = vM(:,7);  bin1656 = vM(:,8);  bin1657 = vM(:,9);
bin1658 = vM(:,10); bin1659 = vM(:,11); bin1660 = vM(:,12); bin1661 = vM(:,13);
bin1662 = vM(:,14); bin1663 = vM(:,15); bin1664 = vM(:,16); bin1665 = vM(:,17);

% Calculate energy
ESD = bin1650.*bin1650 + bin1651.*bin1651 + bin1659.*bin1659 + bin1660.*bin1660...
     +bin1652.*bin1652 + bin1653.*bin1653 + bin1657.*bin1657 + bin1661.*bin1661;
%      +bin1650.*bin1650 + bin1651.*bin1651 + bin1659.*bin1651 + bin1660.*bin1660
%      +bin1650.*bin1650 + bin1651.*bin1651 + bin1659.*bin1651 + bin1660.*bin1660;

Fs = length(ESD)/xtic(length(ESD));
PlotStart = floor(Fs * 37);

% Get rid the impact of BackGround.
BGMK = 10;
UpdataRate = 0;% [0,1]
%calculate model.
BackGroundOld = (1/BGMK)*sum(ESD(1:BGMK));
%update model.
for i = BGMK:1:length(ESD)
    BackGroundNew = (1 - UpdataRate) * BackGroundOld + ESD(i) * UpdataRate;
    UpdataRate = abs(ESD(i) - BackGroundOld) / max(ESD(i), BackGroundOld);
    BackGroundOld = BackGroundNew;
    ESD(i) = BackGroundNew;
end

% p_ESD = ESD(PlotStart:1:length(ESD));
% EsdMean = mean(p_ESD);
% p_ESD= p_ESD(:) - EsdMean;

% Plot the picture.
plot(xtic , ESD );