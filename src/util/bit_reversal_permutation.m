function cols = bit_reversal_permutation(N)
    % 生成位反转排列 (Bit Reversal Permutation)
    % 用于确定极化码的列交换规则
    cols = 0:(N-1);
    n = log2(N);
    for i = 1:N
        x = i - 1;
        val = 0;
        for j = 1:n
            val = bitor(val, bitshift(bitand(x, 1), n - j));
            x = bitshift(x, -1);
        end
        cols(i) = val + 1;
    end
end
