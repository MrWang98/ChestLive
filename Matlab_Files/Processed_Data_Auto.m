close all;
clear;clc;
%%
% 1. This file can be used to automatically execute CA2.m, and help to find the most
% proper arguments. 

% 2. Problematic file nums are recorded in the file - 'problem_nums-vowel_end-name.txt'

% 3. Processed data is stored in the dir named by the object's name.
% Numbers(1, 2, 3...) inside the dir represent the epoch time.(Vowel_end Start from 0.05, increased by 0.02 each epoch.) 

%% 选择数据所在的文件夹
% 设置数据存放的文件夹路径
Path = 'Downloads';      
File = dir(fullfile(Path,'CollectedData_*')); % 显示文件夹下所有符合前缀名为'CollectedData_'的文件夹
all_dir_names = {File.name};
size0 = size(all_dir_names);
length1 = size0(2); % 获取将要处理的文件夹个数
cnt = 1;
fprintf('总文件夹个数： %d \n', length1);

while cnt <= length1
this_dir_name = string(all_dir_names{cnt});
% argument settings
name = this_dir_name; % 处理完毕后的数据所放的文件夹名
data_dir = '/Users/guanqianyun/Downloads/' + name + '/WAV'; % data dir

% 获取音频文件个数
Files = dir(fullfile(data_dir,'SingleWav*'));  
ALLFilesName = {Files.name};
size1 = size(ALLFilesName);
sample_num = size1(2); % sample number
fprintf('Dir %s Start!\n', name);
% sample_num = 20; % sample number

result_root_dir =  '/Users/guanqianyun/Desktop/Processed_Data/'; % result dir
epoch_nums = 15; % 循环的次数
epoch_steps = 0.02 ; % 每次循环增加vowel_end 的大小

% auto-modified arguments
start_vowel = 0.05; % start from 0.05
deleteseg_num = 0; % 0 1 2 3
nxl_list = [ ]; % store the temporal results for debuging
cycle_times = 0; % flag
last_file_path = ''; % 读取上一次有问题的文件编号
first_tag = 1;

%% use 0.05, 0.08, 0.11... process data 
for epoch = 1 : epoch_nums
    vowel_end = start_vowel + (epoch-1) * epoch_steps; % set vowel length
    file_path = sprintf('%s%s%s%s%d%s%s%s', result_root_dir, name, '/', 'problem_nums-', vowel_end, '-', name, '.txt'); % set record file
    file_array_index = 1; % index
    
    if first_tag == 1 % 第一次循环处理
       all_times_number = (file_array_index:1:sample_num); % 文件编号数组
       last_file_path = file_path;
       first_tag = 0;
    
    else % 不是第一次, 则读取上一次有问题文件
        try
            data = textread(last_file_path, '%d');
            all_times_number = data; % 否则记录将要处理的文件标号
        catch
            disp('CONGRATUATIONS, ALL FILES HAVE BEEN SUCCESSFULLY PROCESSED!') %$ 如果不存在文件，说明已经没有异常, 退出程序
            break
        end
        last_file_path = file_path;
    end
    
    fprintf('%s%d  \n', 'Epoch: ', epoch);
    fprintf('%s%d  \n', 'all file number: ', length(all_times_number));
    pause(1)
    
    while file_array_index <= length(all_times_number)
        try
            % Load wav file.
            times_number = all_times_number(file_array_index);
            path_voice1 = '/SingleWav';
            path_voice2 = '.wav';
            path_voice3 = sprintf('%s%s', data_dir, path_voice1);
            path_voice = sprintf('%s%d%s',path_voice3, times_number, path_voice2);
            Data_number = times_number;
            % Read File
            [Rt, Fs] = audioread(path_voice);
             % 获取存储dir，如果不存在则创建
            % store_dir = sprintf('%s%s%s%d%s',result_root_dir, name, '/', epoch, '/');
            store_dir = sprintf('%s%s%s%d%s',result_root_dir, name, '/');
            if exist(store_dir,'dir') == 0
                mkdir(store_dir)
            end
            
            handle_arr = [
                0 % 声音的原始信号。1
                0 % 滤波之后的信号。2
                0 % 希尔伯特之后的包络。3
                0 % 能量谱。4
                1 % MFcc。5
                1 % 音基周期检测。6
                1 % 处理起伏数据。7
                1 % xcorr 8 
                1 % crosscor 9
                1 % Store Data. 10
                0 % Store Voice 11
                0 % reverse 12 default is 0
                0 % Store ESD_pca to .csv 13
            ];

            handle_esd = [
                1 % 分帧 求时间
                1 % BackGround Remove
                1 % EWT
            ];

            %%
            RtBack = Rt;

            %% Delete the begining and end phase.
            RtBack = RtBack(4*Fs:end); 
            Rt = Rt(4*Fs:end); 
            %% Extract audible voice.
            time = length(Rt)/Fs;
            xtic = 0 : 1/Fs :time - 1/Fs;
            [b,a]=butter(1,[200/Fs 1200/Fs],'bandpass');
            Rt=filter(b,a,Rt) * 100;
            Rt=Rt-mean(Rt);    % 消去直流分量
            Rt=Rt/max(abs(Rt));  % 幅值归一化
            % Draw signals
            if (handle_arr(1) == 1)
                figure(1);
                plot(xtic, Rt);
                title("声音的可听见信号");
            end
            %Frequency Features
            N=length(Rt); %计算原信号的长度
            f=Fs*(0:N-1)/N;  %频率分布
            y=fft(Rt);  %对原时域信号x进行fft，得到频域信号y
            if (handle_arr(4) == 1)
                Energy = abs(y) .* abs(y);
                figure(2)
                plot(xtic,Energy)   %绘制滤波之后的时域信号x1
                title('信号的能量谱')
            end
            %设计滤波器FIR滤波器
            b=fir1(8,[200/Fs 1200/Fs]);  %设计带通滤波器 48
            c=freqz(b,1,N);   %频率特性
            y1 = y.*c;
            x1=ifft(y1);
            if (handle_arr(2) == 1)
                 figure(3);
                 plot(xtic,real(x1))   %绘制滤波之后的时域信号x1
                 title('滤波之后的时域信号')
            end
            %hilbert变换，对x1求包络线
            x2=hilbert(real(x1));  %x1的希尔伯特变换x2 real(x1)
            x3=abs(x2);      %x2取模，得到x3
            x4 = - abs(hilbert(real(x1)));
            if (handle_arr(3) == 1)
                figure(4)
                plot(xtic,x3, xtic,x4);
                title('使用希尔伯特之后求得得包络')
            end

            %%
            % 计算音基长度
            if (handle_arr(6) == 1)
                wlen=200;   % 帧长
                inc=80;     % 帧移
                T1= vowel_end;    % 设置基音端点检测的参数 0.05 动态参数 0.5 1
                [voiceseg,vosl,SF,Ef]=pitch_vad(Rt,wlen,inc,T1);   % 基音的端点检测
                fn=length(SF);
                time1 = (0 : length(Rt)-1)/Fs;                % 计算时间坐标
                frameTime = FrameTimeC(fn, wlen, inc, Fs);  % 计算各帧对应的时间坐标
                %figure(5);
                %plot(time1,Rt,'k');  hold on;
                %title('语音信号提取'); axis([0 max(time1) -1 1]);
                %ylabel('幅值'); xlabel('时间/s');
                Result_voice = 0;

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
                %%%%%
                for k=1 : vosl    % 标出有话段
                    nx1=voiceseg(k).begin;
                    nx2=voiceseg(k).end;
                    nxl=voiceseg(k).duration; % 原始的语音长度
                    nxl_list = [nxl_list, nxl]; % 将结果存入数组，为寻找报错原因
                    voice_tmp = Rt(frameTime(nx1) * Fs: frameTime(nx2) * Fs);
                    if (k==1)
                        Result_voice = voice_tmp; % 使用方式Result_voice{1,1}
                    else
                        Result_voice = {Result_voice voice_tmp};
                    end
                    fprintf('%4d   %4d   %4d   %4d\n',k,frameTime(nx1),frameTime(nx2),nxl);
            %         line([frameTime(nx1) frameTime(nx1)],[-1 1],'color','r','linestyle','-');
            %         line([frameTime(nx2) frameTime(nx2)],[-1 1],'color','b','linestyle','--');
                end

                % 将提取的语言绘制出来。

            %     offline = 0;
            %     figure(6)
            %     for i = 1:vosl  % 增加了441的间隔
            %         r_len = length(Result_voice{1,i});
            %         plot((0:441)/Fs + offline,zeros(442),'k');hold on;
            %         r_time = (0 : r_len-1)/Fs + offline + 441/Fs;
            %         plot(r_time, Result_voice{1,i});hold on;
            %         offline = max(r_time);
            %     end
            %     title('绘制语音段');
            %     ylabel('幅值'); xlabel('时间/s');
            end

            %% 
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
                %figure(105) ;
                %ax105 = gca;
                %plot(ax105, frameTime,ccc2(:,1:cn));  % 画出每通道的MFCC系数  1:cn/2
            %     title('(b)MFCC系数');
            %     ylabel('幅值'); xlabel(['时间/s' ]);
            %     figure(8); % 标出有话段部分的 MFCC 。 ================================================================
                ccc2xue = cell(1,2);
                for k=1 : vosl    % 标出有话段
                    nx1=voiceseg(k).begin;
                    nx2=voiceseg(k).end;
                    %line(ax105,[frameTime(nx1) frameTime(nx1)],[-10 10],'color','r','linestyle','-');
                    %line(ax105, [frameTime(nx2) frameTime(nx2)],[-10 10],'color','b','linestyle','--');
                    %ax17_k = subplot(1,2,k);
                    for i = 1:39    % 1:39
                        tmpccc2xue = ccc2(:,i);           
                        if i==1   %获取有胸腔起伏段，语音信号的MFCC
                            ccc2xue{1,k} = tmpccc2xue(nx1:nx2);
                        else
                            ccc2xue{1,k} = [ccc2xue{1,k} tmpccc2xue(nx1:nx2)];
                        end
                        %plot(ax17_k, frameTime(nx1:nx2),tmpccc2xue(nx1:nx2));hold on;
                    end
                end
            end

            %%
            if (handle_arr(7) == 1)
                LenData = length(RtBack); % the length of wav
                time = (0 : LenData-1) / Fs; 
            %    figure(9);plot(time, RtBack); title('Acoustic');
                shape = size(RtBack);
                fftRealArray = zeros(shape)'; 
                LenBatch = 3584; % android microphone buffer size is 3584
                Amp = zeros(ceil(LenData/LenBatch), 4096); %amp contains all the bins 2048

                for n = 1:(LenData/LenBatch)
                    [z,p,k] = butter(9,18000/(Fs/2),'high');
                    [sos_filter ,g] = zp2sos(z,p,k);
                    fftRealArray = RtBack( ((n-1)*3584 + 1) : (n*3584) ); % Extract Batch data.
                    fftRealArray = sosfilt(sos_filter,fftRealArray) *g;
                    % Regulate in amplitude
                    for i = 1:LenBatch
                        fftRealArray(i) = fftRealArray(i) / 32768; 
                    end
                    %Haming Window
                    for i = 1:(LenBatch/2) 
                    %Calculate & apply window symmetrically around center point
                    %Hanning (raised cosine) window
                        winval = 0.5 + 0.5 * cos( (pi*i) / (LenBatch/2));
                        if i > LenBatch/2  
                            winval = 0;
                        end
                        fftRealArray(LenBatch/2 + i) = fftRealArray(LenBatch/2 + i) * winval;
                        fftRealArray(LenBatch/2 + 1 - i) = fftRealArray(LenBatch/2 + 1 - i) * winval;
                    end

                    %fft
                    result = fft(fftRealArray, 4096*2); %4096*2 
                    amp = abs(result);
                    amp = amp(1:4096); 
                    Amp(n,:) = amp' .* 1000;
                end
                [Amp_r, Amp_c] = size(Amp);
                EFs = (Amp_r -1) / ((LenData-1)/Fs);
                ETic = 0:1/EFs:(Amp_r-1)/EFs;
                f_bin = (0:Amp_r-1)*EFs / Amp_r;
            %     figure(10);plot(ETic,Amp(:, 3715));title(['Bins:', num2str(3715)]);xlabel('time (s)'); %1858 

                SubEsd = zeros(22,Amp_r);
                j = 1; time_Features = [];
                for k=3705:1:3726 %3705:1:3725 1850:1:1865
            %         bin_Y = fft(Amp(:, i));
                    bin_Y = smooth(Amp(:, k),10)';
                    if (handle_esd(2) == 1)
                        Esd = bin_Y;
                        %Esd = smooth(Esd,10);
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
            %             figure(k);plot(ETic , Esd );title("Denoise Breathing");%13
                    end
                    SubEsd(j,:) = Esd;
                    j = j+1;

                end
                EsdData = SubEsd(1:j-1,:); %Extract Bins
                % Calcute Cosine Distance.
                y_distance_tmp = pdist(EsdData,'cosine');
                y_distance = squareform(y_distance_tmp); 
                % Calcute idex in all cosine distance.
                y_distance = 1-y_distance;
                [y_distance_sum, y_distance_pos] = max(sum(y_distance));
                y_distance_futher = y_distance(:,y_distance_pos);
                % Remove the fake bins.
                [y_distance_row, y_distance_col] = find(y_distance_futher<0.96); % row is the index.
                rows2remove=unique(y_distance_row); 
                y_distance(rows2remove,:)=[]; 
                % use Kpca
                no_dims =1;
                ESD = 0;
                [mappedD, mapping]= compute_mapping(squeeze(EsdData)', 'ProbPCA',no_dims); 
                KPCA = mappedD;
            %     [mappedD]= compute_mapping(mappedD, 'KPCA',no_dims); 
            %     if (max(mappedD)>=0 && min(mappedD)<=0)
            %         KPCA = -(mappedD-max(mappedD));
            %     else
            %         KPCA = mappedD;  
            %     end
                if(handle_arr(12) == 1)
                    KPCA = -(mappedD-max(mappedD));
                end
                % use Esd
                for i = 1:1:j-1
                    ESD = ESD + SubEsd(i,:) .* SubEsd(i,:);
                end
                FlattenedData_ESD = ESD(:)';
                MappedFlattened_ESD = mapminmax(FlattenedData_ESD, 0, 1);
                MappedData_ESD = reshape(MappedFlattened_ESD, size(ESD)); 
                FlattenedData_KPCA = KPCA(:)';
                MappedFlattened_KPCA = mapminmax(FlattenedData_KPCA, 0, 1);
                MappedData_KPCA = reshape(MappedFlattened_KPCA, size(KPCA)); 
                % 合并
                % data_g = data_gan(times_number, :)';
                % KPCA_front = KPCA(1 : size(KPCA)-size(data_g));
                % size(KPCA_front)
                % KPCA = [KPCA_front', data_g']';
            %     figure(11);plot(ETic, MappedData_ESD);title('ESD');xlabel('time (s)'); 
            %     figure(12);plot(ETic, MappedData_KPCA);title('PPCA ');xlabel('time (s)');
                %% Backgroud remove (Empirical wavelet transform)


                %% EWT Decomposition (Empirical wavelet transform)
                if (handle_esd(3) == 1)
                    %% Initial " signal", "axis", "sampling rate".
                    Esd = KPCA;
                    M = movmean(Esd, 10*EFs); 
                    Esd = Esd - M;
                    FlattenedDataEwt = Esd(:)';
                    MappedFlattenedEwt = mapminmax(FlattenedDataEwt, 0, 1);
                    MappedDataEwt = reshape(MappedFlattenedEwt, size(Esd)); 

                    signal = MappedDataEwt;
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
                    params.InitBounds = [2 25]; % [2 25]
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
            %         figure(14444);
            %         plot(ETic, rec, 'LineWidth',2); 
                end

                %% Locate relate movement
                if( handle_esd(1) == 1)
                    wlen=200;   % 帧长
                    inc=80;     % 帧移
                    y=enframe(Esd,wlen,inc)';    % 分帧
                    fn=size(y,2);                % 求帧数
                    Result_undulate = 0;
                    for k=1 : vosl               % 标出有话段
                        nx1=voiceseg(k).begin;
                        nx2=voiceseg(k).end;
                        Undulate_tmp = rec(floor(frameTime(nx1) * EFs): ceil(frameTime(nx2) * EFs)); %ceil
                        if (k==1)
                            Result_undulate = Undulate_tmp; % 使用方式Result_voice{1,1}
                        else
                            Result_undulate = {Result_undulate Undulate_tmp};
                        end
            %             line([frameTime(nx1) frameTime(nx1)],[-1 1],'color','r','linestyle','-');
            %             line([frameTime(nx2) frameTime(nx2)],[-1 1],'color','b','linestyle','--');
                    end

                end

                %% Extract relate movement
            %     if( handle_esd(1) == 1)
            %         offline_undulate = 0;
            %         figure(15)
            %         for i = 1:vosl  % 增加了441的间隔
            %             r_len_un = length(Result_undulate{1,i});
            %             r_time_un = (0 : r_len_un-1)/EFs + offline_undulate;
            %             plot(r_time_un, Result_undulate{1,i});hold on;
            %             offline_undulate = max(r_time_un);
            %         end
            %         title('绘制起伏段');
            %         ylabel('幅值'); xlabel('时间/s');
            %     end
            end

            %% 计算 xcor
            if (handle_arr(8) == 1)
                %corrcoef xcorr crosscorr
                voiceCundualte1 = xcorr(Result_voice{1,1},Result_undulate{1,1},'none'); % Original cross-correlation calculation
                voiceCundualte2 = xcorr(Result_voice{1,2},Result_undulate{1,2},'none');
            %     voiceCundualte1 = xcorr(voiceGQY.voiceGQY,Result_undulate{1,1},'none'); % UseToMimic Hi
            %     voiceCundualte2 = xcorr(voiceGQY.voiceGQYSiri,Result_undulate{1,2},'none'); % UseToMimic

            %     figure(14);
            %     subplot(1,2,1)
            %     plot(voiceCundualte1);
            %     subplot(1,2,2)
            %     plot(voiceCundualte2);
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
            %     figure(15) ;
            %     subplot(121);
            %     plot(frameTime1,cor1_ccc2(:,1:cor_cn1))  % 画出每通道的MFCC系数  1:cor_cn1  =================1
            %     title('(a)MFCC系数');
            %     ylabel('幅值'); xlabel(['时间/s' ]);
            %     subplot(122);
            %     plot(frameTime2,cor2_ccc2(:,1:cor_cn2))  % 画出每通道的MFCC系数  1:cor_cn2  ==================1
            %     title('(b)MFCC系数');
            %     ylabel('幅值'); xlabel(['时间/s']);
                %===========
                HeySiri = [cor1_ccc2; cor2_ccc2];
            end

            %%
            % 计算 crosscor
            if (handle_arr(9) == 1)
                [acor1,lag1] = crosscorr(Result_voice{1,1},Result_undulate{1,1});
                [acor2,lag2] = crosscorr(Result_voice{1,2},Result_undulate{1,2});
                [~,I1] = max(abs(acor1));
                [~,I2] = max(abs(acor2));
                lagDiff1 = lag1(I1);
                timeDiff1 = lagDiff1/EFs;
                lagDiff2 = lag2(I2);
                timeDiff2 = lagDiff2/EFs;

            %     figure(16)
            %     subplot(211)
            %     plot(lag1,acor1)
            %     a3 = gca;
            %     a3.XTick = sort([-3000:1000:3000 lagDiff1]);
            % 
            %     subplot(212)
            %     plot(lag2,acor2)
            %     a4 = gca;
            %     a4.XTick = sort([-3000:1000:3000 lagDiff2]);
            end 

            %%
            % Store Data
            if (handle_arr(10) == 1)
                % 构建存储文件
                path_store1 = 'HeySiri';
                path_store2 = '.mat';
                path_store_chest = sprintf('%s%s%d%s', store_dir, path_store1, Data_number, path_store2);
                % 文件存在则删除原文件
                if exist(path_store_chest,'file') == 1
                    delete(path_store_chest)
                end
                save(path_store_chest, 'cor1_ccc2',  'cor2_ccc2');
            end

            if (handle_arr(11) == 1)
                voice_store1 = Result_voice{1,1};
                voice_store2 = Result_voice{1,2};
                voicePath_store0  = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/ZhuXiaoTian/Matlab/VoiceStore';
                voicePath_store2 = '.mat';
                path_store_voice = sprintf('%s%d%s',voicePath_store0,Data_number,voicePath_store2);
                save(path_store_voice, 'VoiceStore1', 'VoiceStore2')
                save('voiceStore.mat', 'VoiceStore1', 'VoiceStore2')
            end
            
            if (handle_arr(13) == 1)
                % Store ESD-PCA data
                filename0 = '/Users/guanqianyun/Desktop/Processed_Data/GZ/ESD_pca';
                filename2 = '.csv';
                filename = sprintf('%s%d%s',filename0, Data_number,filename2);
                fid = fopen(filename, 'w');
                for j=1:length(MappedData_KPCA)
                    fprintf(fid, ['%d','\n'], MappedData_KPCA(j));
                end    
            end  
        catch error% 执行出现错误，则调整deleteseg_num
            disp(error)
            % 出现机器无法改正的错误，则输出fail,人工查看错误
            if cycle_times == 1 || length(nxl_list) == 1 
                % print fail
                fprintf('File %4d failed.\n\n',times_number);
                % record the num
                fp = fopen(file_path,'a');
                fprintf(fp, '%d\n', times_number);
                fclose(fp);
                % skip to the next file
                file_array_index = file_array_index + 1;
                cycle_times = 0;
                deleteseg_num = 0;
                nxl_list = [];
                continue;
            end
            % 寻找无效数据
            min_num = inf;
            for i = 1:length(nxl_list)
                if nxl_list(i) < 75 && nxl_list(i) < min_num % 找到最小数据删除
                    deleteseg_num = i;
                    min_num = nxl_list(i);
                end
            end
            % 全部为有效数据,且有三段，则删除最大值
            if deleteseg_num == 0 && length(nxl_list) == 3
                    [value, deleteseg_num] = max(nxl_list); 
            end
            cycle_times = 1; % set flag
            continue;
        end

        fprintf('File %4d has been successfully processed.\n\n',times_number);
        if deleteseg_num ~= 0
            nxl_list(deleteseg_num)= []; % delete the segment number    
        end
        if min(nxl_list) < 75 % if the rest values are also invalid, record the file num.
            % record the num
            fp = fopen(file_path,'a');
            fprintf(fp, '%d\n', times_number);
            fclose(fp);
        end
        % skip to the next file.
        file_array_index = file_array_index + 1;
        cycle_times = 0;
        deleteseg_num = 0;
        nxl_list = [];
    end  
end

cnt = cnt + 1;
end
