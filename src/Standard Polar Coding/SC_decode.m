function decoded_message = SC_decode(codeword, N, K, alpha)
    % 输入:
    % codeword: 编码后的码字 (1xN 的行向量)
    % N: 码长
    % K: 信息比特数量
    % alpha: BEC信道的删除率
    % 输出:
    % decoded_message: 解码后的信息比特 (1xK 的行向量)

    % 初始化解码结果
    decoded_message = zeros(1, K);
    LLR_tree = zeros(log2(N) + 1, N);  % LLR 的解码树

    % 初始 LLR 使用接收的码字
    LLR_tree(end, :) = calculate_initial_llr(codeword, alpha);

    % 递归解码树
    decoded_bits = zeros(1, N);
    decoded_message = decode_node(LLR_tree, decoded_bits, 1, N, 1, K);
end