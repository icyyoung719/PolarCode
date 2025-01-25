function llr = calculate_llr(L, G, i, N)
    % 在 BEC 信道中，LLR 计算应考虑删除的比特
    if L(i) == -1  % 假设 -1 表示该比特被删除
        llr = Inf;  % 删除比特的 LLR 可以是无穷大
    else
        % 对于未删除的比特，计算 LLR
        llr = 2 * sum(L .* G(i, :)) - 0.5;
    end
end
