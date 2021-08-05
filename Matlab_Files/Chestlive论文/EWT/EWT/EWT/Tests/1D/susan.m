%load data
% load S1_GPi_only_data
% D = struct2cell(Data);

%window size
%samp_freq = D{1};
samp_freq =100;
freq_res = 1;
N = samp_freq/freq_res;

%windowing signal
% Trial = D{2};
% info = cell2mat(Trial(1));
% f = info(1,:);

%finding boundaries
V = {};
left = 0;
right = N;
count = 0;
while right < (length(f))
    count = count+1;
    four = fft(f(left+1:right));
    spect = abs(four(1:(length(four)/2)));
    V{count,1} = GSS_BoundariesDetect(spect,'otsu');
    left = left + N/2;
    right = right + N/2;
end

M = zeros(length(V),N/2);
Bound = zeros(length(V),N/2);
for jj = 1:length(V)
    Bound(jj,V{jj}) = V{jj};
end

M(Bound~=0)=1;
B=Bound(Bound~=0);

figure(1);imshow(M,[])
figure(2);histogram(B(:));