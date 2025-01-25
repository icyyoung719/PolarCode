function L = update_L(L, G, i, decoded_bit, N)
    % 更新 L 向量
    % 假设 decoded_bit 是 0 或 1
    if decoded_bit == 1
        L(i) = -L(i);  % 反向 LLR
    end
    % 根据极化码的特性更新 L 向量
    % 具体的更新规则会依赖于极化码的生成矩阵 G
end
