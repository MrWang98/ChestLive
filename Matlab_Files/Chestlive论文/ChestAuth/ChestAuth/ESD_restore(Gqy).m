[x, Fs] = audioread('getWav.wav'); % ??wav??
fprintf('frequecy: %d\n', Fs); % ????

% X: ??; Y: wav???
subplot(2,2,1);
N = length(x);
time = (0 : N-1) / Fs; %???????????
plot(time, x);
xlabel('time');
title('Raw Wav Signal');


% ?????
syms len len_seg fftRealArray Amp
len = length(x); % wav??????
len_seg = 3584; % ??????????????????????segment??
fftRealArray = zeros(3584, 1)'; % ???????segment??fft????
Amp = zeros(len/len_seg, 2048); % ?????fft??????????

% ???segment?????
for n = 1:(len/len_seg)
    [b, a] = butter(1, [19000/Fs 21000/Fs], 'bandpass');
    fftRealArray = x( ((n-1)*3584 + 1) : (n*3584) );
    fftRealArray = filter(b, a, fftRealArray);
    
    % ??????
    for i = 1:len_seg
        fftRealArray(i) = fftRealArray(i) / 32768; 
    end
    
    %Haming???
    for i = 1:(len_seg/2) 
        %Calculate & apply window symmetrically around center point
        %Hanning (raised cosine) window
        winval = 0.5 + 0.5 * cos( (pi*i) / (len_seg/2) );
        if i > len_seg/2  
            winval = 0;
        end
        fftRealArray(len_seg/2 + i) = fftRealArray(len_seg/2 + i) * winval;
        fftRealArray(len_seg/2 + 1 - i) = fftRealArray(len_seg/2 + 1 - i) * winval;
    end
    
    %??FFT?????????Amp?
    result = fft(fftRealArray, 4096);
    amp = abs(result);
    amp = amp(1:2048); 
    Amp(n,:) = amp' .* 1000000;
end

%c?????segment ?ESD?
bin1856 = Amp(:, 1856);bin1857 = Amp(:, 1857);bin1858 = Amp(:, 1858);
bin1859 = Amp(:, 1859);bin1860 = Amp(:, 1860);bin1861 = Amp(:, 1861);
ESD = bin1856 .* bin1856 +  bin1857 .* bin1857 +   bin1858 .* bin1858 + bin1859 .* bin1859 +   bin1860 .* bin1860 + bin1861 .* bin1861  ;

%Movmean get the moving average of an array
M = movmean(ESD, 89);
subplot(2,2,2);
hold on
plot(ESD);
plot(M, 'r');
xlabel('ESD');
ylabel('AMP');
title('Raw ESD');
hold off


% smooth ESD
ESD = smooth(ESD,10);
M = movmean(ESD, 89);
subplot(2,2,3);
hold on
plot(ESD);
plot(M,'r');
xlabel('ESD');
ylabel('AMP');
title('Smooth ESD ')
hold off

% Get rid the impact of BackGround
% To normalize ESD, I modified the process of denoising
BGMK = 10;
UpdataRate = 0;% [0,1]
%calculate model.
BackGround = (1/BGMK) * sum(ESD(1:BGMK));
%old ESD
ESD_old = 0;
%update model.
for i = BGMK:1:length(ESD)
    UpdataRate = abs(ESD(i) - ESD_old) / max(ESD(i), ESD_old); % update rate
    BackGround = (1 - UpdataRate) * BackGround + ESD(i) * UpdataRate; %update model
    ESD_old = ESD(i);
    ESD(i) = ESD(i) - BackGround; % calculate ESD /subtract the Background
end

%draw ESD
M = movmean(ESD, 89);
subplot(2,2,4);
hold on
plot(ESD);
plot(M,'r');
xlabel('ESD');
ylabel('AMP');
title('Smooth ESD with normalized average');
hold off


