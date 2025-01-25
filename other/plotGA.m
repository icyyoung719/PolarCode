% 高斯构造的信道图
N=512;
sigma=0.5; % 表示高斯白噪声的标准差

u = zeros(1, N);
% 初始化高斯信道质量
u(1) = 2/sigma^2;
% 高斯近似递归计算
for i = 1:log2(N)
    j = 2^(i - 1);
    for k = 1:j
        tmp = u(k);
        u(k) = phi_inverse(1 - (1 - phi(tmp))^2);
        u(k + j) = 2 * tmp;
    end
end
u = bitrevorder(u);

scatter((1:N),u(1:N),'.b');
axis([0 1.1*N 0 4*N]);
xlabel('Channel index');
ylabel('E(LLRi)');

% 信道可靠性从左到右逐渐分化，分裂出可靠的信道和不可靠的信道。
% 越高的点表示信道越可靠，可以作为信息位。
% 越低的点表示信道越差，应作为冻结位。
% 信道极化的体现：部分信道非常可靠，部分信道非常差，中间过渡较少。
