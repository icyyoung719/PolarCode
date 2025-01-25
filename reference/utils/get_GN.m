function FN = get_GN(N)
%   N - 极化码的码长，必须是 2 的幂次方（N = 2^n）
%   FN - 生成矩阵 G_N，大小为 N x N

% G_N = F^{⊗n}，其中 F 是基础矩阵

F = [1, 0; 
     1, 1];

FN = zeros(N, N);

% 初始化生成矩阵的第一步，设置 FN 的前 2x2 子矩阵为 F
FN(1 : 2, 1 : 2) = F;

% 通过迭代计算 G_N（逐步进行 Kronecker 乘积扩展）
for i = 2 : log2(N)
    % 取当前矩阵 FN 的前 2^(i-1) x 2^(i-1) 子矩阵
    % 与基础矩阵 F 进行 Kronecker 乘积，更新为 2^i x 2^i 子矩阵
    FN(1 : 2^i, 1 : 2^i) = kron(FN(1 : 2^(i - 1), 1 : 2^(i - 1)), F);
end
end
