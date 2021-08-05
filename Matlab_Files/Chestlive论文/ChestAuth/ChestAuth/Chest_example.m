clear; close all;
% s = settings;
% s.matlab.fonts.codefont.Name.TemporaryValue = 'Times New Roman';
User1 = load('ZhuXiaoTian22_example1.mat'); 
User2 = load('OuRunMin22_example1.mat'); 
User3 = load('QinDang10_example1.mat'); 
U1_Fs = User1.Fs_undulate;
U2_Fs = User2.Fs_undulate;
U3_Fs = User3.Fs_undulate;
U1_selected = User1.Esd(23*U1_Fs:43*U1_Fs-1);
FlattenedData1 = U1_selected(:)';
MappedFlattened1 = mapminmax(FlattenedData1, 0, 1);
MappedData1 = reshape(MappedFlattened1, size(U1_selected)); 
U2_selected = User2.Esd(20*U2_Fs:40*U2_Fs-1);
FlattenedData2 = U2_selected(:)';
MappedFlattened2 = mapminmax(FlattenedData2, 0, 1);
MappedData2 = reshape(MappedFlattened2, size(U2_selected)); 
U3_selected = User3.Esd(16*U3_Fs:36*U3_Fs-1);
FlattenedData3 = U3_selected(:)';
MappedFlattened3 = mapminmax(FlattenedData3, 0, 1);
MappedData3 = reshape(MappedFlattened3, size(U3_selected)); 
U1_xtic = User1.ETic(1:20*U1_Fs);
U2_xtic = User2.ETic(1:20*U2_Fs);
U3_xtic = User3.ETic(1:20*U3_Fs);
figure(1)
h1 = subplot(311);
plot(U1_xtic, MappedData1, 'LineWidth',2); %MappedData3 U1_selected
grid on;
lgd1=legend('U1');
set(gca,'xticklabel',[],'linewidth',1.5,'FontName','Times New Roman','fontsize',22,'GridLineStyle','--');
% title('User1');
h2 = subplot(312);
plot(U2_xtic, MappedData2, 'LineWidth',2);
lgd2=legend('U2');
set(gca,'xticklabel',[],'linewidth',1.5,'FontName','Times New Roman','fontsize',22,'GridLineStyle','--');
% title('User2');
grid on;
ylabel("Normalized Amplitude",'fontsize',22,'FontName','Times New Roman');
h3 = subplot(313);
plot(U3_xtic, MappedData3, 'LineWidth',2);
xlabel("Time (s)",'fontsize',22,'FontName','Times New Roman');
grid on;
% title('User3');
lgd3=legend('U3','FontName','Times New Roman');
set(gca,'linewidth',1.5,'fontsize',22,'FontName','Times New Roman','GridLineStyle','--');


