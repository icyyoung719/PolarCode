function decoded_bits = decode_node(LLR_tree, decoded_bits, level, N, pos, K)
    % 解码树节点的递归过程
    % level: 当前节点所在的树层
    % pos: 当前节点在该层的位置

    if level == log2(N) + 1
        % 叶子节点，直接使用接收码字
        decoded_bits(pos) = 0; % 假设冻结位为 0
    else
        % 递归解码左右子节点
        left_pos = 2 * pos - 1;
        right_pos = 2 * pos;

        % 计算左右子节点的 LLR
        L_left = (1 + LLR_tree(level + 1, left_pos) * LLR_tree(level + 1, right_pos)) / ...
                 (LLR_tree(level + 1, left_pos) + LLR_tree(level + 1, right_pos));
        L_right = L_left;  % 对称性保持简化

        % 更新树的 LLR
        LLR_tree(level, pos) = L_left;

        % 解码左子树
        decode_node(LLR_tree, decoded_bits, level + 1, N, left_pos, K);

        % 利用解码位更新右子树 LLR
        decoded_bit = decoded_bits(left_pos);
        LLR_tree(level + 1, right_pos) = LLR_tree(level + 1, right_pos) * (1 - 2 * decoded_bit);

        % 解码右子树
        decode_node(LLR_tree, decoded_bits, level + 1, N, right_pos, K);
    end
end