close all;
clear;clc;
%加载wav格式的录音文件 /Users/mengxue/Downloads/getWav(1).wav
path = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/TheFirstData/Xusixiao/20200402Vowel/20200402-Vowel-CUP-2.wav';
[Rt, Fs] = audioread(path);

% sound(Rt,Fs);
handle_arr = [
    0 % 声音的原始信号。1
    1 % 滤波之后的信号。2
    0 % 希尔伯特之后的包络。3
    0 % 能量谱。4
    1 % MFcc。5
    1 % 音基周期检测。6
    ];

%添加一个带通滤波器 20000 pilot 确保是高频的声音信号。
time = length(Rt)/Fs;
xtic = 0 : 1/Fs :time - 1/Fs;
[b,a]=butter(1,[200/Fs 1200/Fs],'bandpass');
Rt=filter(b,a,Rt) * 100;
Rt=Rt-mean(Rt);                                % 消去直流分量
Rt=Rt/max(abs(Rt));                            % 幅值归一化

if (handle_arr(1) == 1)
    figure(101);
    plot(xtic, Rt);
    title("声音的原始信号");
end

%Frequency Features
N=length(Rt); %计算原信号的长度
f=Fs*(0:N-1)/N;  %频率分布
y=fft(Rt);  %对原时域信号x进行fft，得到频域信号y
if (handle_arr(4) == 1)
    Energy = abs(y) .* abs(y);
    figure(104)
    plot(xtic,Energy)   %绘制滤波之后的时域信号x1
    title('信号的能量谱')
end

%设计滤波器FIR滤波器
b=fir1(48,[200/Fs 1200/Fs]);  %设计带通滤波器
c=freqz(b,1,N);   %频率特性
y1 = y.*c;
x1=ifft(y1);
if (handle_arr(2) == 1)
    figure(102);
    plot(xtic,real(x1))   %绘制滤波之后的时域信号x1
    title('滤波之后的时域信号')
end

%hilbert变换，对x1求包络线
x2=hilbert(real(x1));  %x1的希尔伯特变换x2 real(x1)
x3=abs(x2);      %x2取模，得到x3
x4 = - abs(hilbert(real(x1)));
if (handle_arr(3) == 1)
    figure(103)
    plot(xtic,x3, xtic,x4);
    title('使用希尔伯特之后求得得包络')
end

% 计算梅尔倒谱系数
if (handle_arr(5) == 1)
    wlen=200; % 帧长
    inc=80;  % 帧移
    num=8; %Mel滤波器个数
    AA = Rt; %Rt real(x1)
    AA=AA/max(abs(AA));    % 幅值归一化 
    time=(0:length(AA)-1)/Fs;
    ccc1=Nmfcc(AA,Fs,num,wlen,inc);
    fn=size(ccc1,1)+4;  %前后各有两帧被舍弃
    cn=size(ccc1,2);
    z=zeros(1,cn);
    ccc2=[z;z;ccc1;z;z];
    frameTime=FrameTimeC(fn,wlen,inc,Fs);               % 求出每帧对应的时间
    figure(105) ; 
    plot(frameTime,ccc2(:,1:cn/2))  % 画出每通道的MFCC系数  1:cn/2
    title('(b)MFCC系数');
    ylabel('幅值'); xlabel(['时间/s' ]);  
end

% 计算音基长度
if (handle_arr(6) == 1)
    wlen=200;   % 帧长
    inc=80;     % 帧移
    T1=0.1;    % 设置基音端点检测的参数 0.05 动态参数 0.5 1
    [voiceseg,vosl,SF,Ef]=pitch_vad(Rt,wlen,inc,T1);   % 基音的端点检测
    fn=length(SF);
    time1 = (0 : length(Rt)-1)/Fs;                % 计算时间坐标
    frameTime = FrameTimeC(fn, wlen, inc, Fs);  % 计算各帧对应的时间坐标
    figure(106);
    plot(time1,Rt,'k');  hold on;
    title('语音信号提取'); axis([0 max(time1) -1 1]);
    ylabel('幅值'); xlabel('时间/s');
    Result_voice = 0;
    for k=1 : vosl    % 标出有话段
        nx1=voiceseg(k).begin;
        nx2=voiceseg(k).end;
        nxl=voiceseg(k).duration; % 原始的语音长度
        voice_tmp = Rt(frameTime(nx1) * Fs: frameTime(nx2) * Fs);
        if (k==1)
            Result_voice = voice_tmp; % 使用方式Result_voice{1,1}
        else
            Result_voice = {Result_voice voice_tmp};
        end
        fprintf('%4d   %4d   %4d   %4d\n',k,frameTime(nx1),frameTime(nx2),nxl);
        line([frameTime(nx1) frameTime(nx1)],[-1 1],'color','r','linestyle','-');
        line([frameTime(nx2) frameTime(nx2)],[-1 1],'color','b','linestyle','--');
    end
    
    % 将提取的语言绘制出来。
    offline = 0;
    figure(107)
    for i = 1:vosl  % 增加了441的间隔
        r_len = length(Result_voice{1,i});
        plot((0:441)/Fs + offline,zeros(442),'k');hold on;
        r_time = (0 : r_len-1)/Fs + offline + 441/Fs;
        plot(r_time, Result_voice{1,i});hold on;
        offline = max(r_time);
    end
    title('绘制语音段');
    ylabel('幅值'); xlabel('时间/s');
end
