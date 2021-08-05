% ZhuXiaoTian   who has the best performance. 
% KuangRuiLing 
% WongJiaSong
% BiHongLiang
% OuRunMin
% DengYangTao 1，5not good.
% LiLingWei
% WuYuan
% QingDang
% TEST_BINS3.csv' ChestLive4.csv'; can get NormalizedEstimation;
% TEST_BINS10 NoPerson
% . ChestInterferenceReduce use WuYuan 74
close all, clear; % StyleTest 
% path_undulate2 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Final/WuYuan/CollectedData/BINS文件/TEST_BINS';
path_undulate2 = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Person/TEST_BINS';
path_undulate3 = '.csv';
times_number = 10; 
path_undulate = sprintf('%s%d%s',path_undulate2,times_number,path_undulate3);
% path_undulate = '/Users/mengxue/Downloads/TEST_BINS10.csv';
BinsSet = csvread(path_undulate, 2,0);

EsdProcess = [
    1 % 计算能量
    1 % 滤除背景噪声。1
    0 % 进行EMD分解。2 We not use EMD.
    0 % 进行EWT
    ];

time = BinsSet(:,1);
xtic= (time(:) - time(1))./1000;
Fs_undulate = length(time)/xtic(length(time)); 
Esd = 0;



% %绘制 胸腔起伏的曲线图
% SubEsd = zeros(21,length(time));
% j = 1; 
% for i= 1:1:16 %1850:1:1865
%     bin_Y = fft(BinsSet(:, i));
%     [bin_x, bin_yt] = max( abs(bin_Y(5:end, 1)) );
%     bin_y = (4+bin_yt -1)*Fs_undulate / length(time);
%     if( bin_y > 0 && bin_y < 10 ) %0.16-0.6
%         SubEsd(j,:) =BinsSet(:, i);
%         j = j+1;
%     end
% end
% EsdData = SubEsd(1:j-1,:); %Extract Bins
% no_dims =1;
% mappedD = Esd
% [Esd, mapping]= compute_mapping(squeeze(EsdData)', 'KPCA',no_dims); 
for i = 8:1:12 % 1856-1860
    Esd = Esd + BinsSet(:,i) .* BinsSet(:,i);
end
Esd_Data_Raw = Esd(36*Fs_undulate:56*Fs_undulate-1);
figure(111112)
plot(xtic,Esd);

figure(12)
MA_Raw = movmean(Esd, 10*Fs_undulate); 
MA_Raw = MA_Raw(36*Fs_undulate:56*Fs_undulate-1);
Esd_static = Esd(36*Fs_undulate:56*Fs_undulate-1);
xtic_static = xtic(1:20*Fs_undulate);
plot(xtic_static,Esd_static);hold  on;
plot(xtic_static,MA_Raw,'r');
title('胸腔起伏的曲线图');
axis([0 max(xtic) -1 1]); 
hold off

% 
% 
% smooth ESD
Esd = smooth(Esd,10);
% Get rid the impact of BackGround
% To normalize ESD, I modified the process of denoising
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


% 
%draw ESD
M = movmean(Esd, 10*Fs_undulate); %89
Esd = Esd - M;
FlattenedData2 = Esd(:)';
MappedFlattened2 = mapminmax(FlattenedData2, 0, 1);
MappedData2 = reshape(MappedFlattened2, size(Esd)); 
figure(222222);
plot(xtic,Esd);hold  on;
plot(xtic,M,'r');

Esd_afterStatic = Esd(36*Fs_undulate:56*Fs_undulate-1); % for drawing ,can delete.
% Draw static effect remove.
figure(123456)
subplot(211);
plot(xtic_static,Esd_static, 'LineWidth',1);hold  on;
grid on;
plot(xtic_static,MA_Raw,':r', 'LineWidth',2);
legend('Original ESD','Moving average of original ESD','FontName','Times New Roman');
ylabel("Amplitude",'fontsize',22,'FontName','Times New Roman');
set(gca,'linewidth',1.5,'FontName','Times New Roman','fontsize',18,'GridLineStyle','--');

subplot(212);
plot(xtic_static,Esd_afterStatic, 'LineWidth',1);hold  on; %before is Esd_afterStatic Esd
grid on;
M = M(36*Fs_undulate:56*Fs_undulate-1);
plot(xtic_static,M,':r', 'LineWidth',2);
legend('ESD after static cancellation','Moving average of original ESD after static cancellation','FontName','Times New Roman');
ylabel("Amplitude",'fontsize',22,'FontName','Times New Roman');
xlabel("Time (s)",'fontsize',22,'FontName','Times New Roman');
set(gca,'linewidth',1.5,'FontName','Times New Roman','fontsize',18,'GridLineStyle','--');

% End Draw

Esd = MappedData2;
Esd = Esd(36*Fs_undulate:56*Fs_undulate-1);
xtic = xtic(1:20*Fs_undulate);
% figure(2222221)
% plot(xtic,Esd_afterStatic);hold  on; %before is Esd_afterStatic Esd
% M = M(36*Fs_undulate:56*Fs_undulate-1);
% plot(xtic,M,'r');
% MappedData2 = Esd; %This line must comment out for ChestBreathNoPerson.m

if (EsdProcess(3) == 1)
    KPCA = MappedData2;
    EFs = Fs_undulate;
    ETic = xtic;
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
    KPCA = MappedData2;
    ETic = xtic;
    EFs = Fs_undulate;
    signal = KPCA;
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
            Show_EWT(ewt,Esd_Data_Raw,rec,EWT_Tic); % add new parameters before Esd_Data_Raw = KPCA
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
%             Show_EWT(ewt,signal,rec);
%             figure(1222); plot(EWT_Tic,signal); % plot the signal before remove other component.
%             title('Before EWT');
%             figure(1333); plot(EWT_Tic,rec); % plot the reconstructed signal.
%             title('After EWT');
        else
             Show_EWT(ewt);
        end    
    end

    if TFplane==1 %Show the time-frequency plane by using the Hilbert transform
        EWT_TF_Plan(ewt,boundaries,params.SamplingRate,signal,[],[],subresf,[]);
    end
    MappedData2 = rec;
    
end

% 
figure(1333)
st_no = 200; et_no = 220; % use to draw No breath. st_no = 1 + 10; et_no = 50 -18; 200; et_no = 300;
esd_start_No = 0:1/Fs_undulate:(et_no*Fs_undulate - st_no*Fs_undulate -1)/Fs_undulate;
esd_end_No = MappedData2(st_no*Fs_undulate:et_no*Fs_undulate -1);
plot(esd_start_No,esd_end_No , 'LineWidth',2); % rec  Esd
xlabel("Time (s)",'fontsize',22,'FontName','Times New Roman');
ylabel("Amplitude",'fontsize',22,'FontName','Times New Roman');
grid on;
set(gca,'linewidth',1,'FontName','Times New Roman','fontsize',22,'GridLineStyle','--');

plot(xtic,MappedData2, 'LineWidth',2); % rec  Esd
hold on

% plot(xtic,M,'r');
% xlabel('Time');
% ylabel('Esd');
% legend('ChestMove','Move Average');
% title('Smooth ESD with normalized average');
% hold off

