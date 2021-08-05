clear; close all;
vad = load('/Users/mengxue/Downloads/TmpData/voice_Vad.mat');
xtic = vad.xtic;
voice = vad.sourcewav;
voiceend = vad.sourceend;
result = vad.vad;
Fs = 16000;
Vad_xtic = xtic(1:20*Fs);
voice_seg = voice(16*Fs:36*Fs-1); % 8 -28
result_seg = result(16*Fs:36*Fs-1);
plot(Vad_xtic, voice_seg);hold on; %voice(1:voiceend)
plot(Vad_xtic, result_seg,'LineWidth',2);
xlabel("Time (s)",'fontsize',22);
ylabel("Normalized ampitude",'fontsize',22,'FontName','Times New Roman');
grid on;
legend('Arouse words','Prediction','fontsize',22,'FontName','Times New Roman');
set(gca,'linewidth',1.5,'FontName','Times New Roman','fontsize',22,'GridLineStyle','--');