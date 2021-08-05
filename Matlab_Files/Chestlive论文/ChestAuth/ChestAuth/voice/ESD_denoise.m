clear,clc;
close all;

% Read data from files (1)
vM = csvread('/Users/mengxue/Documents/Paper/ChestAuthentication/Material/Attackmimic/WithBreath/XueMeng/20200511HiSiri/TEST_BINS2.csv', 2,0);
time = vM(:,1);
xtic= (time(:) - time(1))./1000;
bin1650 = vM(:,2);  bin1651 = vM(:,3);  bin1652 = vM(:,4);  bin1653 = vM(:,5);
bin1654 = vM(:,6);  bin1655 = vM(:,7);  bin1656 = vM(:,8);  bin1657 = vM(:,9);
bin1658 = vM(:,10); bin1659 = vM(:,11); bin1660 = vM(:,12); bin1661 = vM(:,13);
bin1662 = vM(:,14); bin1663 = vM(:,15); bin1664 = vM(:,16); bin1665 = vM(:,17);


% Calculate energy
% Code for calculate ESD
%     for i = 1:1:j-1
%         ESD = ESD + SubEsd(i,:) .* SubEsd(i,:);
%     end
%     figure(13);plot(ETic, ESD);title('ESD');xlabel('time (s)');
ESD = bin1650.*bin1650 + bin1651.*bin1651 + bin1652.*bin1652 + bin1653.*bin1653...
    + bin1654.*bin1654 + bin1657.*bin1657 + bin1655.*bin1655 + bin1656.*bin1656;
figure;plot(xtic , ESD );title("RAW Breathing");

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
figure;plot(xtic , ESD );title("Denoise Breathing");

[imf,residual,info] = emd(ESD);
% emd_visu(ESD,xtic,imf)
f = (0:length(ESD)-1)*Fs/length(ESD);
Imf_freq = 0;
BranchValue = 0;
for i = [2,3,4]
    Y = fft(imf(i,:));
    [x,y] = max(abs(Y));
    if(f(y) > 0.16 && f(y) < 0.6)
        Imf_freq =  f(y);
        BranchValue = i;
    end
end
figure;plot(xtic ,imf(BranchValue,:));title(['EMD branch Value: ', num2str(Imf_freq)]);
