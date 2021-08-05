%MFCC计算函数
function ccc=Nmfcc(x,fs,p,frameSize,inc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 function ccc=Nmfcc(x);
% x是输入语音序列，Mel滤波器的个数为p，采样频率为fs，frameSize为帧长和FFT点数，inc为帧移；ccc为MFCC参数。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 按帧长为frameSize，Mel滤波器的个数为p，采样频率为fs
% 提取Mel滤波器参数，用汉明窗函数
bank=melbankm(p,frameSize,fs,0,0.5,'m');
% 归一化Mel滤波器组系数
bank=full(bank);
bank=bank/max(bank(:));

Mfcc_levels = 13; %阶数的特征

% DCT系数,12*p 
for k=1:Mfcc_levels
  n=0:p-1;
  dctcoef(k,:)=cos((2*n+1)*k*pi/(2*p));
end

% 归一化倒谱提升窗口
w = 1 + 6 * sin(pi * [1:Mfcc_levels] ./ Mfcc_levels);
w = w/max(w);

% 预加重滤波器
xx=double(x);
xx=filter([1 -0.9375],1,xx);

% 语音信号分帧
xx=enframe(xx,frameSize,inc);
n2=fix(frameSize/2)+1;
% 计算每帧的MFCC参数
for i=1:size(xx,1)
  y = xx(i,:);
  s = y' .* hamming(frameSize);
  t = abs(fft(s));
  t = t.^2;
  c1=dctcoef * log(bank * t(1:n2));
  c2 = c1.*w';
  m(i,:)=c2';
end

%差分系数
dtm = zeros(size(m));
for i=3:size(m,1)-2
  dtm(i,:) = -2*m(i-2,:) - m(i-1,:) + m(i+1,:) + 2*m(i+2,:);
end
dtm = dtm / 3;

%求取二阶差分系数
dtmm=zeros(size(dtm));
for i=3:size(dtm,1)-2
dtmm(i,:)=-2*dtm(i-2,:)-dtm(i-1,:)+dtm(i+1,:)+2*dtm(i+2,:);
end
dtmm=dtmm/3;

%合并MFCC参数和一阶差分MFCC参数
ccc = [m dtm dtmm];
%去除首尾两帧，因为这两帧的一阶差分参数为0
ccc = ccc(3:size(m,1)-2,:);

