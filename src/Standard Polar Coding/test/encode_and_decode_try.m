N = 8;  % 码长
K = 4;  % 信息比特数量

% message = [1, 0];  % 原始信息比特
message = randi([0, 1], 1, K);  % 生成包含0或1的随机K位比特序列
alpha = 0.5; % BEC信道的删除率

codeword = standard_encode(message, N, K, alpha);

% decoded_message = SC_decode(codeword, N, K, alpha);

% 输出结果
disp('Message:');
disp(message);
disp('Encoded Codeword:');
disp(codeword);
% disp('Decoded Message:');
% disp(decoded_message);
