% clear;
close all;
% clear;
% vital_path = '/Users/mengxue/Files/DetectDrunk/BigScaleExp/Wangxu/WIFI/wangxu7.csv'

vital_path = '/Users/mengxue/Documents/Paper/ChestAuthentication/Material/BigData/Person/ChestLive4.csv';
breamp = csvread(vital_path, 7,1);

FlattenedData1 = breamp(:)';
MappedFlattened1 = mapminmax(FlattenedData1, 0, 1);
MappedData1 = reshape(MappedFlattened1, size(breamp)) *1.5 -0.0725; 

len = length(breamp);
Fs = 5;
all_time = len / Fs;
% new add
DataDif = 36.6; %16
all_time = (len - DataDif * Fs) / Fs;
MappedData1 = MappedData1(DataDif * Fs:end-1);

% xtic = 0 :1/Fs:all_time - 1/Fs;
xtic_nullog = 0 :1/Fs:all_time -1/Fs;

figure(1);
st = 265;et = 305;
% st = 415;et = 470;
% Plot the figure of nullog
% log_start = xtic_nullog(st*Fs:et*Fs);
log_start = 0:1/Fs:(et*Fs - st*Fs -1)/Fs;
log_end = MappedData1(st*Fs:et*Fs -1);
plot(log_start,log_end,'--','LineWidth',2);
hold on;

% esd_start = xtic(st*Fs_undulate:et*Fs_undulate);
esd_start = 0:1/Fs_undulate:(et*Fs_undulate - st*Fs_undulate -1)/Fs_undulate;
esd_end = MappedData2(st*Fs_undulate:et*Fs_undulate -1);
plot(esd_start, esd_end, 'LineWidth',2);
xlabel("Time (s)",'fontsize',22);
ylabel("Nomarlized estimation",'fontsize',22);
legend('Belt-based estimation','ESD-based estimation','fontsize',22);
grid on;
set(gca,'linewidth',1,'FontName','Times New Roman','fontsize',24,'GridLineStyle','--');

%%
% fft_lenY = et*Fs - st*Fs;
% wY = (0:(et*Fs - st*Fs -1)) * Fs/fft_lenY;
% Ybreath1=fft(MappedData1(st*Fs : et*Fs-1));
% fft_lenY2 = et*Fs_undulate - st*Fs_undulate;
% wY2 = (0:(et*Fs_undulate - st*Fs_undulate -1)) * Fs_undulate/fft_lenY2;
% Ybreath2=fft(MappedData2(st*Fs_undulate : et*Fs_undulate-1));
% figure(155);
% Ybreath1(1) = 0;Ybreath2(1) = 0;
% plot(wY,abs(Ybreath1),'--','LineWidth',2);hold on;
% plot(wY2,abs(Ybreath2),'LineWidth',2);
% xlabel("Frequency (Hz)",'fontsize',22,'FontName','Times New Roman');
% ylabel("Amplitude",'fontsize',22,'FontName','Times New Roman');
% legend('Belt-based frequency','ESD-based frequency','fontsize',24,'FontName','Times New Roman');
% grid on;
% set(gca,'linewidth',1,'FontName','Times New Roman','fontsize',24,'GridLineStyle','--');

%%
% figure(2);
% plot(xtic_nullog,MappedData1);


%%
% % Calculate frequency for noperson.
fft_lenNo = et_no*Fs_undulate - st_no*Fs_undulate;
NoPersonFFT=fft(esd_end_No);
NoPersonFFT(1) = 0;
wNo = (0:(et_no*Fs_undulate - st_no*Fs_undulate -1)) * Fs_undulate/fft_lenNo;
figure(3);
plot(wNo,abs(NoPersonFFT),'LineWidth',2);
xlabel("Frequency (Hz)",'fontsize',22,'FontName','Times New Roman');
ylabel("Amplitude",'fontsize',22,'FontName','Times New Roman');
grid on;
set(gca,'linewidth',1,'FontName','Times New Roman','fontsize',24,'GridLineStyle','--');

