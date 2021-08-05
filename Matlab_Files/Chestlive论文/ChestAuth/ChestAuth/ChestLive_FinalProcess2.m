%% ʹ�ýű�ǰ���Ƚ�line11-13��·��������ȷ��
% 1. Ĭ�ϲ��洢���ݣ��洢���ݵĿ���ѡ���� handle_arr�����У�����10��Store data����Ϊ1���ɴ����ݴ洢��
% 2. ����������д�ӡ�����ݣ����һ��������80��˵����ȡ�����������������line9��deleteseg_num����Ϊ��С�еı�š�
close all;
clear ;clc;
%% Auguments setting19��20��60��77 
times_number = 15; % �ɼ����ݵı�� 0-100
vowel_end = 0.11; % �л�������ѡ�������������Χ��0.05 - 0.4������ͼ106������
deleteseg_num = 0; % ɾ����������ԶΣ���ѡ�� 0 1 2 3 
Data_number = times_number; % ���洢���ݱ�ţ�Ĭ���ǲɼ����ݵı�š�
handle_esd = [
    1 % ��֡ ��ʱ�� 1 
    1 % BackGround Remove 2
    1 % EWT 3
    1 % smooth esd. 4
    1 % fixed frequency.  5
    1 % first 3 segment 6
];
handle_arr = [
    0 % ������ԭʼ�źš�1
    1 % �˲�֮����źš�2
    0 % ϣ������֮��İ��硣3
    0 % �����ס�4
    1 % MFcc��5
    1 % �������ڼ�⡣6
    1 % ����������ݡ�7
    1 % xcorr 8 
    1 % crosscor 9
    0 % Store data 10
    0 % Store voice 11
    ];
% Bins�ļ���Wav�ļ�·����·���е����һ�����ļ����ĵ�һ���֣��ڶ������Ǳ�ţ����������Ǻ�׺��DifferentDirection/XueMeng/CollectedDataRight
path_voice2 = 'C:\Users\dell\Desktop\CollectedData_站立\CollectedData\WAV\SingleWav';
path_undulate2 = 'C:\Users\dell\Desktop\CollectedData_站立\CollectedData\BINS\TEST_BINS';
path_store0  = 'C:\Users\dell\Desktop\CollectedData_站立\CollectedData\HeySiri';


%%
%����wav��ʽ��¼���ļ�
path_voice3 = '.wav';
path_undulate3 = '.csv';
path_voice = sprintf('%s%d%s',path_voice2, times_number, path_voice3);
path_undulate = sprintf('%s%d%s',path_undulate2,times_number,path_undulate3);
[Rt, Fs] = audioread(path_voice);
Rt = Rt(5*Fs:end); %�������ݳ��ȡ�
BinsSet = csvread(path_undulate, 2,0);
% BinsSet = BinsSet(50:end,:);
% sound(Rt,Fs);
%���һ����ͨ�˲��� 20000 pilot ȷ���Ǹ�Ƶ�������źš�
time = length(Rt)/Fs;
xtic = 0 : 1/Fs :time - 1/Fs;
[b,a]=butter(1,[200/Fs 1200/Fs],'bandpass');
Rt=filter(b,a,Rt) * 100;

Rt=Rt-mean(Rt);    % ��ȥֱ������
Rt=Rt/max(abs(Rt));  % ��ֵ��һ��

if (handle_arr(1) == 1)
    figure(101);
    plot(xtic, Rt);
    title("������ԭʼ�ź�");
end

%Frequency Features
N=length(Rt); %����ԭ�źŵĳ���
f=Fs*(0:N-1)/N;  %Ƶ�ʷֲ�
y=fft(Rt);  %��ԭʱ���ź�x����fft���õ�Ƶ���ź�y
if (handle_arr(4) == 1)
    Energy = abs(y) .* abs(y);
    figure(104)
    plot(xtic,Energy)   %�����˲�֮���ʱ���ź�x1
    title('�źŵ�������')
end

%����˲���FIR�˲���
b=fir1(8,[200/Fs 1200/Fs]);  %��ƴ�ͨ�˲��� 48
c=freqz(b,1,N);   %Ƶ������
y1 = y.*c;
x1=ifft(y1);
if (handle_arr(2) == 1)
    figure(102);
    plot(xtic,real(x1))   %�����˲�֮���ʱ���ź�x1
    title('�˲�֮���ʱ���ź�')
end

%hilbert�任����x1�������
x2=hilbert(real(x1));  %x1��ϣ�����ر任x2 real(x1)
x3=abs(x2);      %x2ȡģ���õ�x3
x4 = - abs(hilbert(real(x1)));
if (handle_arr(3) == 1)
    figure(103)
    plot(xtic,x3, xtic,x4);
    title('ʹ��ϣ������֮����õð���')
end

%========================
% ������������
if (handle_arr(6) == 1)
    wlen=200;   % ֡��
    inc=80;     % ֡��
    T1= vowel_end;    % ���û����˵���Ĳ��� 0.05 ��̬���� 0.5 1
    [voiceseg,vosl,SF,Ef]=pitch_vad(Rt,wlen,inc,T1);   % �����Ķ˵���
    fn=length(SF);
    time1 = (0 : length(Rt)-1)/Fs;                % ����ʱ������
    frameTime = FrameTimeC(fn, wlen, inc, Fs);  % �����֡��Ӧ��ʱ������
    figure(106);
    plot(time1,Rt,'k');  hold on;
    title('�����ź���ȡ'); axis([0 max(time1) -1 1]);
    ylabel('��ֵ'); xlabel('ʱ��/s');
%     %%%%%
    if vosl > 2
        i_vosl_temp = vosl;
        if(handle_esd(6) == 1)
            i_vosl_temp =1;
        end    
        remian_num = 3;
        % ȡ�������������
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
    Result_voice = 0;
    %%%%%
    for k=1 : vosl    % ����л���
        nx1=voiceseg(k).begin;
        nx2=voiceseg(k).end;
        nxl=voiceseg(k).duration; % ԭʼ����������
        voice_tmp = Rt(frameTime(nx1) * Fs: frameTime(nx2) * Fs);
        if (k==1)
            Result_voice = voice_tmp; % ʹ�÷�ʽResult_voice{1,1}
        else
            Result_voice = {Result_voice voice_tmp};
        end
        fprintf('%4d   %4d   %4d   %4d\n',k,frameTime(nx1),frameTime(nx2),nxl);
        line([frameTime(nx1) frameTime(nx1)],[-1 1],'color','r','linestyle','-');
        line([frameTime(nx2) frameTime(nx2)],[-1 1],'color','b','linestyle','--');
    end
    
    % ����ȡ�����Ի��Ƴ�����
    offline = 0;
    figure(107)
    for i = 1:vosl  % ������441�ļ��
        r_len = length(Result_voice{1,i});
        plot((0:441)/Fs + offline,zeros(442),'k');hold on;
        r_time = (0 : r_len-1)/Fs + offline + 441/Fs;
        plot(r_time, Result_voice{1,i});hold on;
        offline = max(r_time);
    end
    title('����������');
    ylabel('��ֵ'); xlabel('ʱ��/s');
end

%=================================================================
% ����÷������ϵ�� MFCC
if (handle_arr(5) == 1)
    wlen=200; % ֡��
    inc=80;  % ֡��
    num=8; %Mel�˲�������
    AA = Rt; %Rt real(x1)
    AA=AA/max(abs(AA)); % ��ֵ��һ�� 
    time=(0:length(AA)-1)/Fs;
    ccc1=Nmfcc(AA,Fs,num,wlen,inc);
    fn=size(ccc1,1)+4;  %ǰ�������֡������ ����֡һ�ײ�ֲ���Ϊ0
    cn=size(ccc1,2);
    z=zeros(1,cn);
    ccc2=[z;z;ccc1;z;z];
    frameTime=FrameTimeC(fn,wlen,inc,Fs); % ���ÿ֡��Ӧ��ʱ��
    figure(105) ;
    ax105 = gca;
    plot(ax105, frameTime,ccc2(:,1:cn));  % ����ÿͨ����MFCCϵ��  1:cn/2
    title('(b)MFCCϵ��');
    ylabel('��ֵ'); xlabel(['ʱ��/s' ]);
    figure(17); % ����л��β��ֵ� MFCC �� ================================================================
    ccc2xue = cell(1,2);
    for k=1 : vosl    % ����л���
        nx1=voiceseg(k).begin;
        nx2=voiceseg(k).end;
        line(ax105,[frameTime(nx1) frameTime(nx1)],[-10 10],'color','r','linestyle','-');
        line(ax105, [frameTime(nx2) frameTime(nx2)],[-10 10],'color','b','linestyle','--');
        ax17_k = subplot(1,2,k);
        for i = 1:39    % 1:39
            tmpccc2xue = ccc2(:,i);           
            if i==1   %��ȡ����ǻ����Σ������źŵ�MFCC
                ccc2xue{1,k} = tmpccc2xue(nx1:nx2);
            else
                ccc2xue{1,k} = [ccc2xue{1,k} tmpccc2xue(nx1:nx2)];
            end
            plot(ax17_k, frameTime(nx1:nx2),tmpccc2xue(nx1:nx2));hold on;
        end
    end
end

%=================================================================
% �����������
if (handle_arr(7) == 1)
    % Value of bins
    time = BinsSet(:,1);
    ETic= (time(:) - time(1))./1000;
    Fs_undulate = length(time)/ETic(length(time)); 
    Esd = 0;
    
    %���� ��ǻ���������ͼ
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
    %���� ��ǻ���������ͼ
    if(handle_esd(5) == 1)
        Esd = 0;
        Esd = BinsSet(:,8)  ; %+ BinsSet(:,9)
%         for i = 8:1:12 % 1856-1860
%             Esd = Esd + BinsSet(:,i) .* BinsSet(:,i);
%         end
        handle_esd(3) = 0;
    end

    figure(12)
    FlattenedData2 = Esd(:)';
    MappedFlattened2 = mapminmax(FlattenedData2, -1, 1);
    MappedData2 = reshape(MappedFlattened2, size(Esd));
    plot(ETic,MappedData2);
    title('��ǻ���������ͼ');
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
        figure(134);plot(ETic , Esd );title("Locate");
%         save('QinDang10_example_forAllpath2.mat', 'ETic', 'Fs_undulate','Esd');
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
%         save('QinDang10_ChestVoice1.mat', 'voiceseg','AA','ETic','rec','Fs_undulate','frameTime');
    end
    
    
    % ��֡����λ˵��ʱ�����
    if( handle_esd(1) == 1)
        wlen=200;   % ֡��
        inc=80;     % ֡��
        y=enframe(Esd,wlen,inc)';    % ��֡
        fn=size(y,2);                % ��֡��
        Result_undulate = 0;
        for k=1 : vosl               % ����л���
            nx1=voiceseg(k).begin;
            nx2=voiceseg(k).end;
            Undulate_tmp = MappedData2(floor(frameTime(nx1) * Fs_undulate): ceil(frameTime(nx2) * Fs_undulate)); %ceil
            if (k==1)
                Result_undulate = Undulate_tmp; % ʹ�÷�ʽResult_voice{1,1}
            else
                Result_undulate = {Result_undulate Undulate_tmp};
            end
            line([frameTime(nx1) frameTime(nx1)],[-1 1],'color','r','linestyle','-');
            line([frameTime(nx2) frameTime(nx2)],[-1 1],'color','b','linestyle','--');
        end
        
    end

    % ����ȡ��������Ƴ�����
    offline_undulate = 0;
    figure(13)
    for i = 1:vosl  % ������441�ļ��
        r_len_un = length(Result_undulate{1,i});
        r_time_un = (0 : r_len_un-1)/Fs_undulate + offline_undulate;
        plot(r_time_un, Result_undulate{1,i});hold on;
        offline_undulate = max(r_time_un);
    end
    title('���������');
    ylabel('��ֵ'); xlabel('ʱ��/s');
end


%=================================================================
% ���� xcor
if (handle_arr(8) == 1)
    %corrcoef xcorr crosscorr
    voiceCundualte1 = xcorr(Result_voice{1,1},Result_undulate{1,1},'none');
    voiceCundualte2 = xcorr(Result_voice{1,2},Result_undulate{1,2},'none');
    
    figure(14);
    subplot(1,2,1)
    plot(voiceCundualte1);
    subplot(1,2,2)
    plot(voiceCundualte2);
    %===========
    cor_wlen=200; % ֡��
    cor_inc=80;  % ֡��
    cor_num=8; %Mel�˲�������
    
    cor_AA1 = voiceCundualte1; 
    cor_AA2 = voiceCundualte2;
    cor_AA1=cor_AA1/max(abs(cor_AA1)); % ��ֵ��һ�� 
    cor_AA2=cor_AA2/max(abs(cor_AA2)); % ��ֵ��һ�� 
    cor_time1=(0:length(cor_AA1)-1)/Fs;
    cor_time2=(0:length(cor_AA2)-1)/Fs;
    cor1_ccc1=Nmfcc(cor_AA1,Fs,num,wlen,inc);
    cor_fn1=size(cor1_ccc1,1)+4;  %ǰ�������֡������
    cor_cn1=size(cor1_ccc1,2);
    cor_z1=zeros(1,cor_cn1);
    cor1_ccc2=[cor_z1;cor_z1;cor1_ccc1;cor_z1;cor_z1];
    
    cor2_ccc1=Nmfcc(cor_AA2,Fs,num,wlen,inc);
    cor_fn2=size(cor2_ccc1,1)+4;  %ǰ�������֡������
    cor_cn2=size(cor2_ccc1,2);
    cor_z2=zeros(1,cor_cn2);
    cor2_ccc2=[cor_z2;cor_z2;cor2_ccc1;cor_z2;cor_z2];
    frameTime1=FrameTimeC(cor_fn1,wlen,inc,Fs); % ���ÿ֡��Ӧ��ʱ��
    frameTime2=FrameTimeC(cor_fn2,wlen,inc,Fs); % ���ÿ֡��Ӧ��ʱ��
    figure(15) ;
    subplot(121);
    plot(frameTime1,cor1_ccc2(:,1:cor_cn1))  % ����ÿͨ����MFCCϵ��  1:cor_cn1  =================1
    title('(a)MFCCϵ��');
    ylabel('��ֵ'); xlabel(['ʱ��/s' ]);
    subplot(122);
    plot(frameTime2,cor2_ccc2(:,1:cor_cn2))  % ����ÿͨ����MFCCϵ��  1:cor_cn2  ==================1
    title('(b)MFCCϵ��');
    ylabel('��ֵ'); xlabel(['ʱ��/s']);
    %===========
    HeySiri = [cor1_ccc2; cor2_ccc2];
end

%=================================================================
% ���� crosscor
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

%%
% Save Data.
if (handle_arr(10) == 1)
    path_store2 = '.mat';
    path_store_chest = sprintf('%s%d%s',path_store0,Data_number,path_store2);
    save(path_store_chest, 'cor1_ccc2',  'cor2_ccc2');
end

if (handle_arr(11) == 1)
    voice_store1 = Result_voice{1,1};
    voice_store2 = Result_voice{1,2};
    voicePath_store0  = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/Attack/Matlab/VoiceStore';
    voicePath_store2 = '.mat';
    path_store_voice = sprintf('%s%d%s',voicePath_store0,Data_number,voicePath_store2);
    save(path_store_voice, 'voice_store1', 'voice_store2')
%     save('voiceStore.mat', 'VoiceStore1', 'VoiceStore2')
end