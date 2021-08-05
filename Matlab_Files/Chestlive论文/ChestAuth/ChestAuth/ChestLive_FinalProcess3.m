%% 使用脚本前，先将line11-13的路径设置正确。
% 1. 默认不存储数据，存储数据的控制选项在 handle_arr变量中，将第10项Store data设置为1即可打开数据存储。
% 2. 检测命令行中打印的数据，最后一列若大于80则说明提取正常。多大于三，则将line9的deleteseg_num设置为最小行的编号。
close all;
clear ;clc;
%% Auguments setting
times_number = 1; % 采集数据的编号 0-100
vowel_end = 0.05; % 有话语音段选择参数：调整范围，0.05 - 0.4。根据图106调整。
deleteseg_num = 0; % 删除多余的语言段：可选项 0 1 2 3 4
deleteseg_num1 = 0; % 删除多余的语言段：可选项 0 1 2 3 4
Data_number = times_number; % 当存储数据编号，默认是采集数据的编号。

segments_voice = [
    0 % 0 表示三段，1表示4段
    0 
];

handle_esd = [
    1 % 分帧 求时间 1 
    1 % BackGround Remove 2
    1 % EWT 3
    0 % smooth esd. 4
    1 % fixed frequency.  5
    0 % first 3 segment 6
];
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
    0 % Store data 10
    0 % Store voice 11
    ];
% Bins文件和Wav文件路径，路径中的最后一项是文件名的第一部分，第二部分是编号，第三部分是后缀。
path_voice2 = 'C:\Users\dell\Desktop\CollectedData_走路\CollectedData\WAV\SingleWav';
path_undulate2 = 'C:\Users\dell\Desktop\CollectedData_走路\CollectedData\BINS\TEST_BINS';
path_store0  = 'C:\Users\dell\Desktop\CollectedData_走路\CollectedData\HeySiri';


%%
%加载wav格式的录音文件
path_voice3 = '.wav';
path_undulate3 = '.csv';
path_voice = sprintf('%s%d%s',path_voice2, times_number, path_voice3);
path_undulate = sprintf('%s%d%s',path_undulate2,times_number,path_undulate3);
[Rt, Fs] = audioread(path_voice);
Rt = Rt(5*Fs:end); %限制数据长度。
BinsSet = csvread(path_undulate, 2,0);
BinsSet = BinsSet(50:end,:);
% sound(Rt,Fs);
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
    if vosl > 4
        i_vosl_temp = vosl;
        if(handle_esd(6) == 1)
            i_vosl_temp =1;
        end    
        remian_num = 5;
        % 取最后两个语音。
        while(i_vosl_temp-remian_num>0)
            i_vosl_temp1 = i_vosl_temp-remian_num;
            voiceseg(i_vosl_temp1) = [];
            i_vosl_temp = i_vosl_temp - 1;
            if(handle_esd(6) == 1)
                i_vosl_temp = i_vosl_temp + 1;
            end    
        end
        vosl = remian_num;
    end
    
    if(deleteseg_num ~= 0)
        voiceseg(deleteseg_num) = []; %Used temporal. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vosl = vosl -1;           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    if(deleteseg_num1 ~= 0)
        voiceseg(deleteseg_num1) = []; %Used temporal. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        vosl = vosl -1;           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    Result_voice = cell(1,vosl);
    %%%%%
    for k=1 : vosl    % 标出有话段
        nx1=voiceseg(k).begin;
        nx2=voiceseg(k).end;
        nxl=voiceseg(k).duration; % 原始的语音长度
        voice_tmp = Rt(frameTime(nx1) * Fs: frameTime(nx2) * Fs);
        if k==1   %获取有胸腔起伏段，语音信号的MFCC
            Result_voice{1,k} = voice_tmp;
        else
            Result_voice{1,k} = [Result_voice{1,k} voice_tmp];
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
%         plot((0:441)/Fs + offline,zeros(442),'k');hold on;
        plot((0:44)/Fs + offline,zeros(45),'k');hold on;
        r_time = (0 : r_len-1)/Fs + offline + 44/Fs;
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
    ccc2xue = cell(1,vosl);
    for k=1 : vosl    % 标出有话段
        nx1=voiceseg(k).begin;
        nx2=voiceseg(k).end;
        line(ax105,[frameTime(nx1) frameTime(nx1)],[-10 10],'color','r','linestyle','-');
        line(ax105, [frameTime(nx2) frameTime(nx2)],[-10 10],'color','b','linestyle','--');
        ax17_k = subplot(1,vosl,k); % sub figure number.
        for i = 1:39    % 1:39
            tmpccc2xue = ccc2(:,i);           
            if i==1   %获取有胸腔起伏段，语音信号的MFCC
                ccc2xue{1,k} = tmpccc2xue(nx1:nx2);
            else
                ccc2xue{1,k} = [ccc2xue{1,k} tmpccc2xue(nx1:nx2)];
            end
            plot(ax17_k, frameTime(nx1:nx2), tmpccc2xue(nx1:nx2));hold on;
        end
    end
end

%=================================================================
% 处理起伏数据
if (handle_arr(7) == 1)
    % Value of bins
    time = BinsSet(:,1);
    ETic= (time(:) - time(1))./1000;
    Fs_undulate = length(time)/ETic(length(time)); 
    Esd = 0;
    
    %绘制 胸腔起伏的曲线图
    SubEsd = zeros(21,length(time));
    j = 1; 
    for i= 1:1:16 %1850:1:1865
        bin_Y = fft(BinsSet(:, i));
        [bin_x, bin_yt] = max( abs(bin_Y(5:end, 1)) );
        bin_y = (4+bin_yt -1)*Fs_undulate / length(time);
        if( bin_y > 0.16 && bin_y < 0.6 ) % bin_y > 0.16 && bin_y < 0.6
            SubEsd(j,:) =BinsSet(:, i);
            j = j+1;
        end
    end
    EsdData = SubEsd(1:j-1,:); %Extract Bins
    for i = 1:1:j-1
        Esd = Esd + EsdData(i,:).* EsdData(i,:);
    end
    %绘制 胸腔起伏的曲线图
    if(handle_esd(5) == 1)
        Esd = 0;
        for i = 8:1:12 % 1856-1860
            Esd = Esd + BinsSet(:,i) .* BinsSet(:,i);
        end
        handle_esd(3) = 0;
    end

    figure(12)
    FlattenedData2 = Esd(:)';
    MappedFlattened2 = mapminmax(FlattenedData2, -1, 1);
    MappedData2 = reshape(MappedFlattened2, size(Esd));
    plot(ETic,MappedData2);
    title('胸腔起伏的曲线图');
    axis([0 max(ETic) -1 1]); 

        %% Backgroud remove (Empirical wavelet transform)
    if (handle_esd(2) == 1)
        if(handle_esd(4) == 1)
            Esd = smooth(Esd,10)';
        end

        BGMK = 10;
        UpdataRate = 0;% [0,1]
        %calculate model.
        BackGround = (1/BGMK) * sum(Esd(1:BGMK));
        %old ESD
        BackGround_old = 0;
        %update model.
        for i = BGMK:1:length(Esd)
            UpdataRate = abs(Esd(i) - BackGround_old) / max(Esd(i), BackGround_old); % update rate
            BackGround = (1 - UpdataRate) * BackGround + Esd(i) * UpdataRate; %update model
            BackGround_old = Esd(i); 
            Esd(i) = Esd(i) - BackGround; % calculate ESD /subtract the Background
        end
        figure(133);plot(ETic , Esd );title("Denoise Breathing");
    end
    
    %% EWT Decomposition (Empirical wavelet transform)
    if (handle_esd(3) == 1)
        %% Initial " signal", "axis", "sampling rate".
        EFs = Fs_undulate;
        M = movmean(Esd, 10*EFs);
        Esd = Esd - M;
        FlattenedDataEwt = Esd(:)';
        MappedFlattenedEwt = mapminmax(FlattenedDataEwt, 0, 1);
        MappedDataEwt = reshape(MappedFlattenedEwt, size(Esd)); 
       
        signal = MappedDataEwt'; % MappedDataEwt
        EWT_Tic = ETic;
        params.SamplingRate = EFs; %put -1 if you don't know the sampling rate

        %%
        % Choose the wanted global trend removal (none,plaw,poly,morpho,tophat)
        params.globtrend = 'none';
        params.degree=10; % degree for the polynomial interpolation
        % Choose the wanted regularization (none,gaussian,avaerage,closing)
        params.reg = 'none';
        params.lengthFilter = 6;
        params.sigmaFilter = 1.5;
        % Choose the wanted detection method (locmax,locmaxmin,ftc, adaptive,adaptivereg,scalespace)
        params.detect = 'scalespace';
        params.typeDetect='otsu'; %for scalespace:otsu,halfnormal,empiricallaw,mean,kmeans
        params.N = 3; % maximum number of bands
        params.completion = 0; % choose if you want to force to have params.N modes
                               % in case the algorithm found less ones (0 or 1)
        params.InitBounds = [2 25];
        % Perform the detection on the log spectrum instead the spectrum
        params.log=0;
        % Choose the results you want to display (Show=1, Not Show=0)
        Bound=0;   % Display the detected boundaries on the spectrum
        Comp=1;    % Display the EWT components
        Rec=1;     % Display the reconstructed signal
        TFplane=0; % Display the time-frequency plane (by using the Hilbert 
                   % transform). You can decrease the frequency resolution by
                   % changing the subresf variable below.
        subresf=1;
        InitBounds = params.InitBounds;

        [ewt,mfb,boundaries]=EWT1D(signal,params);

        %% Show the results
        if Bound==1 %Show the boundaries on the spectrum
            div=1;
            if (strcmp(params.detect,'adaptive')||strcmp(params.detect,'adaptivereg'))
                Show_EWT_Boundaries(abs(fft(KPCA)),boundaries,div,params.SamplingRate,InitBounds);
            else
                Show_EWT_Boundaries(abs(fft(KPCA)),boundaries,div,params.SamplingRate);
            end
        end

        if Comp==1 %Show the EWT components and the reconstructed signal
            if Rec==1
                %compute the reconstruction
                rec=iEWT1D(ewt,mfb);
    %             Show_EWT(ewt,KPCA,rec);
                EWT_Component_Show = 0;%plot the EWT components
                if EWT_Component_Show ==1
                    figure; l=1;
                    if length(ewt)>6
                        lm=6;
                    else
                        lm=length(ewt);
                    end
                    for k=1:length(ewt)
                        hold on; subplot(lm,1,l); plot(EWT_Tic,ewt{k}); %axis off;
                        if mod(k,6) == 0
                            figure;
                            l=1;
                        else
                            l=l+1;
                        end
                    end
                end
%                 figure(133); plot(EWT_Tic,signal); % plot the signal before remove other component.
%                 title('Before EWT');
%                 figure(1333); plot(EWT_Tic,rec); % plot the reconstructed signal.
%                 title('After EWT');
                 Esd = rec;
            else
                 Show_EWT(ewt);
            end    
        end

        if TFplane==1 %Show the time-frequency plane by using the Hilbert transform
            EWT_TF_Plan(ewt,boundaries,params.SamplingRate,signal,[],[],subresf,[]);
        end
        figure(1444);
        plot(ETic, rec, 'LineWidth',2); 
    end
    
    
    % 分帧，定位说话时的起伏
    if( handle_esd(1) == 1)
        wlen=200;   % 帧长
        inc=80;     % 帧移
        y=enframe(Esd,wlen,inc)';    % 分帧
        fn=size(y,2);                % 求帧数
        Result_undulate = cell(1,vosl); % 起伏段
        for k=1 : vosl               % 标出有话段
            nx1=voiceseg(k).begin;
            nx2=voiceseg(k).end;
            Undulate_tmp = MappedData2(floor(frameTime(nx1) * Fs_undulate): ceil(frameTime(nx2) * Fs_undulate)); %ceil
            if k==1   %获取有胸腔起伏段，语音信号的MFCC
                Result_undulate{1,k} = Undulate_tmp;
            else
                Result_undulate{1,k} = [Result_undulate{1,k} Undulate_tmp];
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
    voiceCundualte1 = xcorr(Result_voice{1,1},Result_undulate{1,1},'none');
    voiceCundualte2 = xcorr(Result_voice{1,2},Result_undulate{1,2},'none');
    voiceCundualte3 = xcorr(Result_voice{1,3},Result_undulate{1,3},'none');
    if(segments_voice(1) == 1)
        voiceCundualte4 = xcorr(Result_voice{1,4},Result_undulate{1,4},'none');
    end
    
    figure(14);
    subplot(1,3,1)
    plot(voiceCundualte1);
    subplot(1,3,2)
    plot(voiceCundualte2);
    subplot(1,3,3)
    plot(voiceCundualte3);
    %===========
    cor_wlen=200; % 帧长
    cor_inc=80;  % 帧移
    cor_num=8; %Mel滤波器个数
    
    cor_AA1 = voiceCundualte1; 
    cor_AA2 = voiceCundualte2;
    cor_AA3 = voiceCundualte3;
    if(segments_voice(1) == 1)
        cor_AA4 = voiceCundualte4;
    end
    cor_AA1=cor_AA1/max(abs(cor_AA1)); % 幅值归一化 
    cor_AA2=cor_AA2/max(abs(cor_AA2)); % 幅值归一化 
    cor_AA3=cor_AA3/max(abs(cor_AA3)); % 幅值归一化 
    if(segments_voice(1) == 1)
        cor_AA4=cor_AA4/max(abs(cor_AA4)); % 幅值归一化 
    end
    cor_time1=(0:length(cor_AA1)-1)/Fs;
    cor_time2=(0:length(cor_AA2)-1)/Fs;
    cor_time3=(0:length(cor_AA3)-1)/Fs;
    if(segments_voice(1) == 1)
        cor_time4=(0:length(cor_AA4)-1)/Fs;% 幅值归一化 
    end
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
    
    cor3_ccc1=Nmfcc(cor_AA3,Fs,num,wlen,inc);
    cor_fn3=size(cor3_ccc1,1)+4;  %前后各有两帧被舍弃
    cor_cn3=size(cor3_ccc1,2);
    cor_z3=zeros(1,cor_cn3);
    cor3_ccc2=[cor_z3;cor_z3;cor3_ccc1;cor_z3;cor_z3];
    if(segments_voice(1) == 1)
        cor4_ccc1=Nmfcc(cor_AA4,Fs,num,wlen,inc);
        cor_fn4=size(cor4_ccc1,1)+4;  %前后各有两帧被舍弃
        cor_cn4=size(cor4_ccc1,2);
        cor_z4=zeros(1,cor_cn4);
        cor4_ccc2=[cor_z4;cor_z4;cor4_ccc1;cor_z4;cor_z4];
    end
    frameTime1=FrameTimeC(cor_fn1,wlen,inc,Fs); % 求出每帧对应的时间
    frameTime2=FrameTimeC(cor_fn2,wlen,inc,Fs); % 求出每帧对应的时间
    frameTime3=FrameTimeC(cor_fn3,wlen,inc,Fs); % 求出每帧对应的时间
    if(segments_voice(1) == 1)
        frameTime4=FrameTimeC(cor_fn4,wlen,inc,Fs); % 求出每帧对应的时间
    end
    figure(15) ;
    subplot(1,vosl,1);
    plot(frameTime1,cor1_ccc2(:,1:cor_cn1))  % 画出每通道的MFCC系数  1:cor_cn1  =================1
    title('(a)MFCC系数');
    ylabel('幅值'); xlabel(['时间/s' ]);
    subplot(1,vosl,2);
    plot(frameTime2,cor2_ccc2(:,1:cor_cn2))  % 画出每通道的MFCC系数  1:cor_cn2  ==================1
    title('(b)MFCC系数');
    ylabel('幅值'); xlabel(['时间/s']);
    subplot(1,vosl,3);
    plot(frameTime3,cor3_ccc2(:,1:cor_cn3))  % 画出每通道的MFCC系数  1:cor_cn2  ==================1
    title('(c)MFCC系数');
    ylabel('幅值'); xlabel(['时间/s']);
    if(segments_voice(1) == 1)
        subplot(1,vosl,4);
        plot(frameTime4,cor4_ccc2(:,1:cor_cn4))  % 画出每通道的MFCC系数  1:cor_cn2  ==================1
        title('(d)MFCC系数');
        ylabel('幅值'); xlabel(['时间/s']);
    end
    %===========
    OKGoogle = [cor1_ccc2; cor2_ccc2;cor3_ccc2];
    if(segments_voice(1) == 1)
        OKGoogle = [cor1_ccc2; cor2_ccc2;cor3_ccc2;cor4_ccc2];
    end
end

%=================================================================
% 计算 crosscor
if (handle_arr(9) == 1)
    [acor1,lag1] = crosscorr(Result_voice{1,1},Result_undulate{1,1});
    [acor2,lag2] = crosscorr(Result_voice{1,2},Result_undulate{1,2});
    [acor3,lag3] = crosscorr(Result_voice{1,3},Result_undulate{1,3});
    [~,I1] = max(abs(acor1));
    [~,I2] = max(abs(acor2));
    [~,I3] = max(abs(acor3));
    lagDiff1 = lag1(I1);
    timeDiff1 = lagDiff1/Fs_undulate;
    lagDiff2 = lag2(I2);
    timeDiff2 = lagDiff2/Fs_undulate;
    lagDiff3 = lag3(I3);
    timeDiff3 = lagDiff3/Fs_undulate;

    figure(16)
    subplot(311)
    plot(lag1,acor1)
    a3 = gca;
%     a3.XTick = sort([-3000:1000:3000 lagDiff1]);

    subplot(312)
    plot(lag2,acor2)
    a4 = gca;
%     a4.XTick = sort([-3000:1000:3000 lagDiff2]);

    subplot(313)
    plot(lag3,acor3)
    a5 = gca;
%     a5.XTick = sort([-3000:1000:3000 lagDiff3]);
end 

%%
% Save Data.
if (handle_arr(10) == 1)
    path_store2 = '.mat';
    path_store_chest = sprintf('%s%d%s',path_store0,Data_number,path_store2);
    save(path_store_chest, 'cor1_ccc2',  'cor2_ccc2', 'cor3_ccc2');
    if(segments_voice(1) == 1)
        save(path_store_chest, 'cor1_ccc2',  'cor2_ccc2', 'cor3_ccc2','cor4_ccc2');
    end
end

if (handle_arr(11) == 1)
    voice_store1 = Result_voice{1,1};
    voice_store2 = Result_voice{1,2};
    voice_store3 = Result_voice{1,3};
    voicePath_store0  = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/ZhuXiaoTian/Matlab/VoiceStore';
    voicePath_store2 = '.mat';
    path_store_voice = sprintf('%s%d%s',voicePath_store0,Data_number,voicePath_store2);
    save(path_store_voice, 'voice_store1', 'voice_store2', 'voice_store3')
%     save('voiceStore.mat', 'VoiceStore1', 'VoiceStore2')
end