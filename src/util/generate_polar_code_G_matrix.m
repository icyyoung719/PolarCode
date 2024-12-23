function G = generate_polar_code_G_matrix(n)
    % 参数:
    %   n - 指定递归层数，决定生成矩阵G的大小。
    % 返回值:
    %   G - 极化码生成矩阵，大小为 2^n x 2^n。


    % 基础矩阵 F
    F = [1, 0; 1, 1];

    % 初始化 G' 为 F
    G_prime = F;

    % 通过 Kronecker 乘积递归构造 G'
    for i = 2:n
        G_prime = kron(G_prime, F);
    end

    cols = bit_reversal_permutation(2^n);
    % 交换列得到真正的 G 矩阵
    G = G_prime(:, cols);
end
