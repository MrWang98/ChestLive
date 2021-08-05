close all;
clear ;clc;
% Load other person's voice
voiceGQY = load('voiceGQY_Cup.mat'); 

%加载wav格式的录音文件
% path_voice1 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/Attackmimic/WithOutBreath/XuSiXiao/20200512Hisiri/';
% path_voice2 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/TheFirstData/Xusixiao/20200408Vowel/VOWEL-CUP-';
path_voice2 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/XueMeng/CollectedData/WAV文件/getWav';
path_voice3 = '.wav';
times_number = 83;
Data_number = times_number;
vowel_end =0.17; % line 104f
deleteseg_num = 0;
path_voice = sprintf('%s%d%s',path_voice2, times_number, path_voice3);
[Rt, Fs] = audioread(path_voice);
% Rt1 = Rt;
% Read Bins data from files
% path_undulate1 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/Attackmimic/WithOutBreath/XuSiXiao/20200512Hisiri/';
path_undulate2 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/XueMeng/CollectedData/BINS文件/TEST_BINS';
path_undulate3 = '.csv';
path_undulate = sprintf('%s%d%s',path_undulate2,times_number,path_undulate3);
BinsSet = csvread(path_undulate, 2,0);
handle_esd = [
    1 % 分帧 求时间
    ];

% sound(Rt,Fs);
handle_arr = [
    0 % 声音的原始信号。1
    1 % 滤波之后的信号。2
    0 % 希尔伯特之后的包络。3
    0 % 能量谱。4
    1 % MFcc。5
    1 % 音基周期检测。6
    1 % 处理起伏数据。7
    1 % xcorr 8 
    1 % crosscor 9
    ];

%添加一个带通滤波器 20000 pilot 确保是高频的声音信号。
time = length(Rt)/Fs;
xtic = 0 : 1/Fs :time - 1/Fs;
[b,a]=butter(1,[200/Fs 1200/Fs],'bandpass');
Rt=filter(b,a,Rt) * 100;

Rt=Rt-mean(Rt);    % 消去直流分量
Rt=Rt/max(abs(Rt));  % 幅值归一化

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
b=fir1(8,[200/Fs 1200/Fs]);  %设计带通滤波器 48
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

%========================
% 计算音基长度
if (handle_arr(6) == 1)
    wlen=200;   % 帧长
    inc=80;     % 帧移
    T1= vowel_end;    % 设置基音端点检测的参数 0.05 动态参数 0.5 1
    [voiceseg,vosl,SF,Ef]=pitch_vad(Rt,wlen,inc,T1);   % 基音的端点检测
    fn=length(SF);
    time1 = (0 : length(Rt)-1)/Fs;                % 计算时间坐标
    frameTime = FrameTimeC(fn, wlen, inc, Fs);  % 计算各帧对应的时间坐标
    figure(106);
    plot(time1,Rt,'k');  hold on;
    title('语音信号提取'); axis([0 max(time1) -1 1]);
    ylabel('幅值'); xlabel('时间/s');
%     %%%%%
    if vosl > 2
        i_vosl_temp = vosl;
        remian_num = 3;
        % 取最后两个语音。
        while(i_vosl_temp-remian_num>0)
            i_vosl_temp1 = i_vosl_temp-remian_num;
            voiceseg(i_vosl_temp1) = [];
            i_vosl_temp = i_vosl_temp - 1;
        end
        vosl = remian_num;
    end
    
    if(deleteseg_num ~= 0)
        voiceseg(deleteseg_num) = []; %Used temporal. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vosl = vosl -1;           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    Result_voice = 0;
    %%%%%
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

%=================================================================
% 计算梅尔倒谱系数 MFCC
if (handle_arr(5) == 1)
    wlen=200; % 帧长
    inc=80;  % 帧移
    num=8; %Mel滤波器个数
    AA = Rt; %Rt real(x1)
    AA=AA/max(abs(AA)); % 幅值归一化 
    time=(0:length(AA)-1)/Fs;
    ccc1=Nmfcc(AA,Fs,num,wlen,inc);
    fn=size(ccc1,1)+4;  %前后各有两帧被舍弃 这两帧一阶差分参数为0
    cn=size(ccc1,2);
    z=zeros(1,cn);
    ccc2=[z;z;ccc1;z;z];
    frameTime=FrameTimeC(fn,wlen,inc,Fs); % 求出每帧对应的时间
    figure(105) ;
    ax105 = gca;
    plot(ax105, frameTime,ccc2(:,1:cn));  % 画出每通道的MFCC系数  1:cn/2
    title('(b)MFCC系数');
    ylabel('幅值'); xlabel(['时间/s' ]);
    figure(17); % 标出有话段部分的 MFCC 。 ================================================================
    ccc2xue = cell(1,2);
    for k=1 : vosl    % 标出有话段
        nx1=voiceseg(k).begin;
        nx2=voiceseg(k).end;
        line(ax105,[frameTime(nx1) frameTime(nx1)],[-10 10],'color','r','linestyle','-');
        line(ax105, [frameTime(nx2) frameTime(nx2)],[-10 10],'color','b','linestyle','--');
        ax17_k = subplot(1,2,k);
        for i = 1:39    % 1:39
            tmpccc2xue = ccc2(:,i);           
            if i==1   %获取有胸腔起伏段，语音信号的MFCC
                ccc2xue{1,k} = tmpccc2xue(nx1:nx2);
            else
                ccc2xue{1,k} = [ccc2xue{1,k} tmpccc2xue(nx1:nx2)];
            end
            plot(ax17_k, frameTime(nx1:nx2),tmpccc2xue(nx1:nx2));hold on;
        end
    end
end

%=================================================================
% 处理起伏数据
if (handle_arr(7) == 1)
    % Value of bins
    time = BinsSet(:,1);
    xtic= (time(:) - time(1))./1000;
    Fs_undulate = length(time)/xtic(length(time)); 
    Esd = 0;
    figure(11)
    % Plot the bins.
    for i = 2:1:17   
        plot(xtic,BinsSet(:,i)); hold on;
    end
    title('不同bins的曲线');

    %绘制 胸腔起伏的曲线图
    for i = 8:1:12 % 1856-1860
        Esd = Esd + BinsSet(:,i) .* BinsSet(:,i);
    end
    figure(12)
    FlattenedData2 = Esd(:)';
    MappedFlattened2 = mapminmax(FlattenedData2, -1, 1);
    MappedData2 = reshape(MappedFlattened2, size(Esd));
    plot(xtic,MappedData2);
    title('胸腔起伏的曲线图');
    axis([0 max(xtic) -1 1]); 

    % 分帧，定位说话时的起伏
    if( handle_esd(1) == 1)
        wlen=200;   % 帧长
        inc=80;     % 帧移
        y=enframe(Esd,wlen,inc)';    % 分帧
        fn=size(y,2);                % 求帧数
        Result_undulate = 0;
        for k=1 : vosl               % 标出有话段
            nx1=voiceseg(k).begin;
            nx2=voiceseg(k).end;
            Undulate_tmp = MappedData2(floor(frameTime(nx1) * Fs_undulate): ceil(frameTime(nx2) * Fs_undulate)); %ceil
            if (k==1)
                Result_undulate = Undulate_tmp; % 使用方式Result_voice{1,1}
            else
                Result_undulate = {Result_undulate Undulate_tmp};
            end
            line([frameTime(nx1) frameTime(nx1)],[-1 1],'color','r','linestyle','-');
            line([frameTime(nx2) frameTime(nx2)],[-1 1],'color','b','linestyle','--');
        end
        
    end

    % 将提取的起伏绘制出来。
    offline_undulate = 0;
    figure(13)
    for i = 1:vosl  % 增加了441的间隔
        r_len_un = length(Result_undulate{1,i});
        r_time_un = (0 : r_len_un-1)/Fs_undulate + offline_undulate;
        plot(r_time_un, Result_undulate{1,i});hold on;
        offline_undulate = max(r_time_un);
    end
    title('绘制起伏段');
    ylabel('幅值'); xlabel('时间/s');
end


%=================================================================
% 计算 xcor
if (handle_arr(8) == 1)
    %corrcoef xcorr crosscorr
%     voiceCundualte1 = xcorr(Result_voice{1,1},Result_undulate{1,1},'none');
    voiceCundualte2 = xcorr(Result_voice{1,2},Result_undulate{1,2},'none');
    voiceCundualte1 = xcorr(voiceGQY.voiceGQY,Result_undulate{1,1},'none'); % UseToMimic Hi
%     voiceCundualte2 = xcorr(voiceGQY.voiceGQYSiri,Result_undulate{1,2},'none'); % UseToMimic
    
    figure(14);
    subplot(1,2,1)
    plot(voiceCundualte1);
    subplot(1,2,2)
    plot(voiceCundualte2);
    %===========
    cor_wlen=200; % 帧长
    cor_inc=80;  % 帧移
    cor_num=8; %Mel滤波器个数
    
    cor_AA1 = voiceCundualte1; 
    cor_AA2 = voiceCundualte2;
    cor_AA1=cor_AA1/max(abs(cor_AA1)); % 幅值归一化 
    cor_AA2=cor_AA2/max(abs(cor_AA2)); % 幅值归一化 
    cor_time1=(0:length(cor_AA1)-1)/Fs;
    cor_time2=(0:length(cor_AA2)-1)/Fs;
    cor1_ccc1=Nmfcc(cor_AA1,Fs,num,wlen,inc);
    cor_fn1=size(cor1_ccc1,1)+4;  %前后各有两帧被舍弃
    cor_cn1=size(cor1_ccc1,2);
    cor_z1=zeros(1,cor_cn1);
    cor1_ccc2=[cor_z1;cor_z1;cor1_ccc1;cor_z1;cor_z1];
    
    cor2_ccc1=Nmfcc(cor_AA2,Fs,num,wlen,inc);
    cor_fn2=size(cor2_ccc1,1)+4;  %前后各有两帧被舍弃
    cor_cn2=size(cor2_ccc1,2);
    cor_z2=zeros(1,cor_cn2);
    cor2_ccc2=[cor_z2;cor_z2;cor2_ccc1;cor_z2;cor_z2];
    frameTime1=FrameTimeC(cor_fn1,wlen,inc,Fs); % 求出每帧对应的时间
    frameTime2=FrameTimeC(cor_fn2,wlen,inc,Fs); % 求出每帧对应的时间
    figure(15) ;
    subplot(121);
    plot(frameTime1,cor1_ccc2(:,1:cor_cn1))  % 画出每通道的MFCC系数  1:cor_cn1  =================1
    title('(a)MFCC系数');
    ylabel('幅值'); xlabel(['时间/s' ]);
    subplot(122);
    plot(frameTime2,cor2_ccc2(:,1:cor_cn2))  % 画出每通道的MFCC系数  1:cor_cn2  ==================1
    title('(b)MFCC系数');
    ylabel('幅值'); xlabel(['时间/s']);
    %===========
    HiSiri = [cor1_ccc2; cor2_ccc2];
end

%=================================================================
% 计算 crosscor
if (handle_arr(9) == 1)
    [acor1,lag1] = crosscorr(Result_voice{1,1},Result_undulate{1,1});
    [acor2,lag2] = crosscorr(Result_voice{1,2},Result_undulate{1,2});
    [~,I1] = max(abs(acor1));
    [~,I2] = max(abs(acor2));
    lagDiff1 = lag1(I1);
    timeDiff1 = lagDiff1/Fs_undulate;
    lagDiff2 = lag2(I2);
    timeDiff2 = lagDiff2/Fs_undulate;

    figure(16)
    subplot(211)
    plot(lag1,acor1)
    a3 = gca;
%     a3.XTick = sort([-3000:1000:3000 lagDiff1]);

    subplot(212)
    plot(lag2,acor2)
    a4 = gca;
%     a4.XTick = sort([-3000:1000:3000 lagDiff2]);
end 

%=================================================================
% Save Data.
% path_store0  = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/TheFirstData/Xusixiao/MimicGQY/CUP';
% % path_store1  = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/Attackmimic/WithOutBreath/XuSiXiao/Matlab20200512/HiSiri';
% path_store2 = '.mat';
% path_store_chest = sprintf('%s%d%s',path_store0,Data_number,path_store2);
% % path_store_voice = sprintf('%s%d%s',path_store0,times_number,path_store2);
% save(path_store_chest, 'cor1_ccc2'); %, 'cor2_ccc2', 'HiSiri'
% % save(path_store_voice, 'ccc2xue'); %  'cor1_ccc2', 'cor2_ccc2' two word， cup，luck。
%================================================================
% voiceGQYHi = Result_voice{1,1};
% voiceGQYSiri = Result_voice{1,2};
% save('voiceGQY.mat', 'voiceGQYHi', 'voiceGQYSiri')