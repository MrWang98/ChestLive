close all;
clear;clc;
%加载wav格式的录音文件
path_voice = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/TheFirstData/Xusixiao/20200402Vowel/20200402-Vowel-CAT-1.wav';
[Rt, Fs] = audioread(path_voice);
%添加一个带通滤波器 20000 pilot 确保是高频的声音信号。
[b,a]=butter(1,[19000/Fs 21000/Fs],'bandpass');
Rt=filter(b,a,Rt);
 
%计算录音文件的ESD，sine频率为21KHz，考虑20hz的多普勒，bins为 976 977 978 
N = 2048;
L = length(Rt);
for i = 1 :N: L-N
    if( i == 1)%对第一个2048点进行FFT，并求ESD  (Rt(1:1:N)是指接收信号的2048个点。 930:932
        j = 1;
        y_tmp = fft((Rt(1:1:N)),N) ;
        Rf = abs(y_tmp) .* abs(y_tmp);
        Esd =  sum(Rf(976:978)) / N; %pilot 为 20khz
    else%对录音文件中剩余的2048点进行FFT，并求ESD
        j = j + 1;
        y_tmp = fft( Rt(N*(j-1)+1 :1: N*j), N);
        Rf = abs(y_tmp) .* abs(y_tmp);
        tmp = sum(Rf(976:978))/N;
        Esd = [Esd tmp];
    end
end

%平滑波形
Esd = smoothdata(Esd);

%根据绘制ESD的时域图  2048/48000
len_esd = length(Esd);
t = 0 :N/Fs: (len_esd - 1) * N/Fs;
figure;
plot(t,Esd);
xlabel("时间：单位 秒s");
ylabel("ESD");