function bit = llr_to_bit(llr, bit_idx, u_hat, N)
    % 使用递归方式计算 SC 的 LLR 更新和判决
    % 输入:
    % llr: 当前的对数似然比向量
    % bit_idx: 当前译码比特的索引
    % u_hat: 当前的译码比特向量
    % N: 当前码字长度
    % 输出:
    % bit: 当前译码比特

    if N == 1
        % 基础情况：仅一个比特，直接根据LLR判决
        bit = llr(bit_idx) < 0;
    else
        % 计算中间结果
        half_N = N / 2;
        if bit_idx <= half_N
            % 左分支计算
            bit = llr_to_bit(f_l(llr(1:half_N), llr(half_N+1:end), u_hat(1:half_N)), bit_idx, u_hat, half_N);
        else
            % 右分支计算
            bit = llr_to_bit(f_g(llr(1:half_N), llr(half_N+1:end), u_hat(1:half_N)), bit_idx - half_N, u_hat, half_N);
        end
    end
end