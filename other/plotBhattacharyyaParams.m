clc;
clear all;

% 初始化参数
indx = 1:10;
N = 2.^indx;
C = zeros(max(N), max(N));
C(1, 1) = 1 - 0.5; % 初始化第一个巴氏参数,其alpha=0.5

% 计算巴氏参数
for idx = 1:length(indx)
    N0 = N(idx);
    for j = 1:N0/2 
        C(N0, 2*j-1) = C(N0/2, j)^2;
        C(N0, 2*j) = 2*C(N0/2, j) - C(N0/2, j)^2;
    end
end

figure;


N0 = 128;
subplot(2, 1, 1); % 创建两个子图中的第一个
plot(1:N0, C(N0, 1:N0), 'k*');
axis([1 N0, 0 1]);
xlabel('i');
ylabel(['C_{', num2str(N0), '}^i']);
title(['Bhattacharyya Parameters for N=', num2str(N0)]);


N0 = 1024;
subplot(2, 1, 2); % 创建两个子图中的第二个
plot(1:N0, C(N0, 1:N0), 'k*');
axis([1 N0, 0 1]);
xlabel('i');
ylabel(['C_{', num2str(N0), '}^i']);
title(['Bhattacharyya Parameters for N=', num2str(N0)]);