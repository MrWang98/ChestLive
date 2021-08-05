
data_real = csvread('real_breath_signals490.csv', 0, 0);

data_g = csvread('generated_ESD490.csv', 0, 0);
Mapped_data_g = mapminmax(data_g, 0, 1); % 归一化原数据

data_esd = csvread('input_ESD490.csv', 0, 0);

index = 3;
figure(1)
% plot(data_real(index, :))
hold on;
plot(Mapped_data_g(index, :))
hold on;
plot(data_esd(index, :))
