function decoded_message = SC_decode(received, N, K, alpha)
    % 输入:
    % received: 接收到的比特向量 (1×N 的行向量)
    % N: 码长 (必须是 2 的幂)
    % K: 信息比特的数量 (K < N)
    % alpha: BEC信道的删除率
    % 输出:
    % decoded_message: 解码后的信息比特 (1×K 的行向量)

    % Step 1: 构造 G 矩阵
    G = generate_polar_code_G_matrix(log2(N));
    
    % Step 2: 用巴氏参数法计算 Z 向量, 存放信道错误率
    Z = channel_capacity_vector(N, alpha);
    
    % Step 3: 根据 Z 向量排序，确定冻结比特的位置
    [~, indices] = sort(Z, 'descend');  % 按照信道错误率降序排序
    frozen_indices = indices(1:N-K);  % 冻结比特的位置
    info_indices = indices(N-K+1:end);  % 信息比特的位置
    
    % Step 4: Successive Cancellation (SC) 解码
    decoded_message = zeros(1, K);
    u_vector = zeros(1, N);  % 初始假设所有比特都是 0


    %{
    
    for i = 1:N
        if ismember(i, info_indices)  % 如果是信息比特
            % 通过 SC 解码推断该比特的值
            % 如果接收到的信号中该比特是错误的，改正之
            if received(i) == 0  % 如果收到的是 0
                u_vector(i) = 0;  % 假设为 0
            else  % 如果收到的是 1
                u_vector(i) = 1;  % 假设为 1
            end
        end
    end
    
    %}
    
    % Step 5: 从解码的 u_vector 中提取出信息比特
    decoded_message = u_vector(info_indices);
end
