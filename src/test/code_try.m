function polar_encoding_bhattacharyya(N, K)
    % Polar编码实现，并使用巴氏参数构造可靠子信道

    % N: 码字长度（必须是2的幂）
    % K: 信息位长度
    
    % 生成极化矩阵
    F = [1 0; 1 1]; % 基本极化矩阵
    G = kron(F, F);  % 生成N维极化矩阵
    
    % 初始化信息比特
    u = zeros(1, N);  % u为信息位（初始化为0）
    info_bits = randi([0 1], 1, K);  % 随机生成K个信息位
    
    % 将信息位插入到u的前K个位置
    u(1:K) = info_bits;
    
    % 生成极化编码输出
    x = polar_encode(u, G, N);
    
    % 计算巴氏参数，评估子信道可靠性
    bhattacharyya = bhattacharyya_parameter(G, N);
    
    % 输出编码结果
    disp('极化编码输出：');
    disp(x);
    
    % 输出巴氏参数
    disp('巴氏参数（子信道可靠性）：');
    disp(bhattacharyya);
end

function x = polar_encode(u, G, N)
    % 极化编码函数
    x = u * G;  % 计算编码结果
    x = mod(x, 2);  % 确保结果是二进制
end

function bhattacharyya = bhattacharyya_parameter(G, N)
    % 计算巴氏参数，用于评估每个子信道的可靠性
    bhattacharyya = zeros(1, N);
    for i = 1:N
        % 假设巴氏参数计算公式为：B(i) = P(0->1) + P(1->0)
        % 这只是一个近似方法，具体计算可以更复杂
        bhattacharyya(i) = 2^(-i);  % 假设值，实际计算可以用更精确的模型
    end
end
