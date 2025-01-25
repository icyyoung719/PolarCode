function u_vector = recursive_SC_decode(l_vector, u_vector, frozen_indices, l, r)
    % 递归 SC 译码函数
    % 输入:
    % l_vector: 接收到的码字 (1×N 的行向量)
    % u_vector: 初始 u 向量 (1×N 的行向量)
    % frozen_indices: 冻结比特的位置 (1×N-K 的行向量)
    % l: 当前段的左边界
    % r: 当前段的右边界
    % 输出:
    % u_vector: 更新后的 u 向量

    if r - l + 1 == 1
        % 基本情况: 只有一个比特
        if ismember(l, frozen_indices)
            u_vector(l) = 0;  % 冻结比特填 0
        else
            u_vector(l) = l_vector(l);  % 信息比特直接取接收的值
        end
    else
        % 递归情况: 分治
        mid = floor((l + r) / 2);

        % Step 1: 译码左半部分
        l1 = l_vector(l:mid);
        l2 = l_vector(mid+1:r);
        l_combined_left = mod(l1 + l2, 2);  % 模 2 相加
        u_vector = recursive_SC_decode(l_combined_left, u_vector, frozen_indices, l, mid);

        % Step 2: 译码右半部分
        l_combined_right = mod(l_vector(mid+1:r) + u_vector(l:mid), 2);  % 模 2 相加
        u_vector = recursive_SC_decode(l_combined_right, u_vector, frozen_indices, mid+1, r);
    end
end