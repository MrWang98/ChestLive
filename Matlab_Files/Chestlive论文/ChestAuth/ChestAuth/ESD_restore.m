%% 使用了降低维工具箱： what drtoolbox/techniques or help compute_mapping
% FactorAnalysis  PCA
% PCA is overed by demension reduce tool by relate path.
close all; clear;
% WavPath = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/Attackmimic/WithBreath/XueMeng/20200511HiSiri/SingleWav19.wav';
% WavPath = '/Users/mengxue/Downloads/getWav124.wav';
WavPath = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/ZhuTianLin/CollectedData/WAV文件/getWav47.wav';
[x, Fs] = audioread(WavPath);
EsdProcess = [
    1 % 计算能量
    1 % 滤除背景噪声。1
    1 % 进行EMD分解。2 We not use EMD.
    0 % 进行EWT
    ];

%%==========Look acoustic wave form==========
LenData = length(x); % the length of wav
time = (0 : LenData-1) / Fs; 
figure(11);plot(time, x); title('Acoustic');

shape = size(x);
fftRealArray = zeros(shape)'; 
LenBatch = 3584; % android microphone buffer size is 3584
Amp = zeros(LenData/LenBatch, 4096); %amp contains all the bins 2048

%%
for n = 1:(LenData/LenBatch)
%     [b, a] = butter(6, [19000/Fs 21000/Fs], 'bandpass'); %1
    [z,p,k] = butter(9,18000/(Fs/2),'high');
    [sos_filter ,g] = zp2sos(z,p,k);
%     y2=sosfilt(sos_filter,x);%滤波
%     x=filter(b,a,x);
    fftRealArray = x( ((n-1)*3584 + 1) : (n*3584) ); % Extract Batch data.
%     fftRealArray = filter(b, a, fftRealArray);
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
    result = fft(fftRealArray, 4096*2); 
    amp = abs(result);
    amp = amp(1:4096); 
    Amp(n,:) = amp' .* 1000;
end
[Amp_r, Amp_c] = size(Amp);
EFs = (Amp_r -1) / ((LenData-1)/Fs);
ETic = 0:1/EFs:(Amp_r-1)/EFs;
f_bin = (0:Amp_r-1)*EFs / Amp_r;
figure(12);plot(ETic,Amp(:, 3715));title(['Bins:', num2str(3715)]);xlabel('time (s)'); %1858 
% xueM = Amp(:, 3715);
% xueM1 = normalize(xueM(1:10));
% figure(1222);plot(xueM1);title(['Bins:', num2str(3715)]);xlabel('time (s)'); %1858 
    
%%  Scan the bins, then calculate Dimensionality reduction
if (EsdProcess(1) == 1)
    SubEsd = zeros(21,Amp_r);
    j = 1; 
    for i=3705:1:3725
        bin_Y = fft(Amp(:, i));
        [bin_x, bin_yt] = max( abs(bin_Y(5:end, 1)) );
        bin_y = (4+bin_yt -1)*EFs / Amp_r;
        if( bin_y > 0 && bin_y < 10.6 ) % bin_y > 0.16 && bin_y < 0.6
            SubEsd(j,:) =Amp(:, i);
            j = j+1;
        end
    end
    EsdData = SubEsd(1:j-1,:); %Extract Bins
    no_dims =1;
%     [mappedD, mapping]= compute_mapping(squeeze(EsdData)', 'KPCA',no_dims); 
%     [mappedD]= compute_mapping(mappedD, 'KPCA',no_dims); 
%     if (max(mappedD)>=0 && min(mappedD)<=0)
%         KPCA = -(mappedD-max(mappedD));
%     end
%     Code for calculate ESD
    ESD = 0;
    for i = 1:1:j-1
        ESD = ESD + SubEsd(i,:) .* SubEsd(i,:);
    end
    figure(13);plot(ETic, ESD);title('ESD');xlabel('time (s)');
%     figure(133);plot(ETic, KPCA);title('KPCA ');xlabel('time (s)');
end

% Get rid the impact of BackGround.
if (EsdProcess(2) == 1)
    KPCA = ESD;
    BGMK = 10;
    UpdataRate = 0;% [0,1]
    %calculate model.
    BackGroundOld = (1/BGMK)*sum(KPCA(1:BGMK));
    %update model.
    for i = BGMK:1:length(KPCA)
        BackGroundNew = (1 - UpdataRate) * BackGroundOld + KPCA(i) * UpdataRate;
        UpdataRate = abs(KPCA(i) - BackGroundOld) / max(KPCA(i), BackGroundOld);
        BackGroundOld = BackGroundNew;
        KPCA(i) = BackGroundNew;
    end
    figure(21);plot(ETic , KPCA );title("Denoise Breathing");
end

%% EMD Decomposition
if (EsdProcess(3) == 1)
    [imf,residual,info] = emd(KPCA);
%     emd_visu(ESD,ETic,imf)
    f = (0:length(KPCA)-1)*EFs/length(KPCA);
    Imf_freq = 0;
    BranchValue = 0;
    for i = [2,3,4]
        Y = fft(imf(i,:));
        [x,y] = max(abs(Y));
        if(f(y) > 0.16 && f(y) < 0.6) % 0.16
            Imf_freq =  f(y);
            BranchValue = i;
        end
    end
    figure(22);plot(ETic ,imf(BranchValue,:));title(['EMD branch Value: ', num2str(Imf_freq)]);
end

%% EWT Decomposition (Empirical wavelet transform)
if (EsdProcess(4) == 1)
    %% Initial " signal", "axis", "sampling rate".
    signal = KPCA;
    EWT_Tic = ETic;
    params.SamplingRate = EFs; %put -1 if you don't know the sampling rate
    
    %%
    % Choose the wanted global trend removal (none,plaw,poly,morpho,tophat)
    params.globtrend = 'none';
    params.degree=6; % degree for the polynomial interpolation
    % Choose the wanted regularization (none,gaussian,avaerage,closing)
    params.reg = 'none';
    params.lengthFilter = 10;
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
            figure(); plot(EWT_Tic,signal); % plot the signal before remove other component.
            title('Before EWT');
            figure(); plot(EWT_Tic,rec); % plot the reconstructed signal.
            title('After EWT');
        else
             Show_EWT(ewt);
        end    
    end

    if TFplane==1 %Show the time-frequency plane by using the Hilbert transform
        EWT_TF_Plan(ewt,boundaries,params.SamplingRate,signal,[],[],subresf,[]);
    end
    
end
