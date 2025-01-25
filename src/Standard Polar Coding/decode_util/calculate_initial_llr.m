function initial_llr = calculate_initial_llr(codeword, alpha)
    % 计算初始 LLR 对应 BEC 信道模型
    initial_llr = zeros(1, length(codeword));
    for i = 1:length(codeword)
        if codeword(i) == 1
            initial_llr(i) = log((1 - alpha) / alpha);
        elseif codeword(i) == 0
            initial_llr(i) = -log((1 - alpha) / alpha);
        else
            initial_llr(i) = 0;  % 擦除
        end
    end
end
