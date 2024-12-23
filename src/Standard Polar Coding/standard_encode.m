function codeword = standard_encode(message, N, K, alpha)
    % 输入:
    % message: 原始信息比特向量 (1×K 的行向量)
    % N: 码长 (必须是 2 的幂)
    % K: 信息比特的数量 (K < N)
    % alpha: BEC信道的删除率
    % 输出:
    % codeword: 编码后的码字 (1×N 的行向量)

    % Step 1: 构造 G 矩阵
    G = generate_polar_code_G_matrix(log2(N));
    
    % Step 2: 用巴氏参数法计算 Z 向量,存放信道错误率
    Z = channel_capacity_vector(N, alpha);
    
    % Step 3: 根据 Z 向量排序，确定冻结比特的位置
    [~, indices] = sort(Z, 'descend');  % 按照信道错误率降序排序
    frozen_indices = indices(1:N-K); 
    info_indices = indices(N-K+1:end); 

    % Step 4: 构造 u_vector，填入信息和冻结比特
    u_vector = zeros(1, N);
    u_vector(info_indices) = message;  
    u_vector(frozen_indices) = 0;  % 冻结比特填 0

    % Step 5: 计算 codeword = u_vector * G
    codeword = mod(u_vector * G, 2); 
end
