clear; close all;
User1 = load('QinDang10_example_forAllpath1.mat'); 
User2 = load('QinDang10_example_forAllpath2.mat'); 
U1_Fs = User1.Fs_undulate;
U2_Fs = User2.Fs_undulate;
U1_selected = User1.Esd(16*U1_Fs:36*U1_Fs-1);
FlattenedData1 = U1_selected(:)';
MappedFlattened1 = mapminmax(FlattenedData1, 0, 1);
MappedData1 = reshape(MappedFlattened1, size(U1_selected)); 
U2_selected = User2.Esd(16*U2_Fs:36*U2_Fs-1);
FlattenedData2 = U2_selected(:)';
MappedFlattened2 = mapminmax(FlattenedData2, 0, 1);
MappedData2 = reshape(MappedFlattened2, size(U2_selected)); 
U1_xtic = User1.ETic(1:20*U1_Fs);
U2_xtic = User2.ETic(1:20*U2_Fs);
figure(1)
plot(U1_xtic, MappedData1, ':','LineWidth',2);hold on;
plot(U2_xtic, MappedData2, 'LineWidth',2);
xlabel("Time (s)",'fontsize',22,'FontName','Times New Roman'); 
ylabel("Normalized ampitude",'fontsize',22,'FontName','Times New Roman'); 
grid on;
legend('Fixed range  ','All range','fontsize',20,'FontName','Times New Roman');
set(gca,'linewidth',1.5,'FontName','Times New Roman','fontsize',18,'GridLineStyle','--');
% figure(2)
% plot(User2.ETic, User2.Esd);


