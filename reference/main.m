clear;

% -------------------------------编码的基本参数------------------------------
crc_len = 4;             % CRC校验位长度
n = 12;                   % 极化码的bit位数
N = 2^n;                 % 极化码长度
R = 0.5;                 % 码率
max_runs = 1000;         % 仿真运行次数
msgbit_len = N * R;      % 发送信息比特的长度
K = msgbit_len + crc_len; % 编码后的总位长度
snr = 2;                 % 高斯信道信噪比（分贝）
L = 16;                  % CASCL译码器的列表长度

% ------------------------------预处理---------------------------------
[gen, det] = get_crc_objective(crc_len); % 输出对应的CRC生成器与CRC校验器

lambda_offset = 2.^(0 : log2(N));  % 分段向量，表示每一层的信道分段长度
llr_layer_vec = get_llr_layer(N);  % LLR计算的实际执行层数向量
bit_layer_vec = get_bit_layer(N);  % 比特值返回时实际执行层数向量

% -------------------------生成信息比特并添加CRC校验------------------------
info = randsrc(msgbit_len, 1, [0 1; 0.5 0.5]); % 等概率生成信息比特
info_with_crc = gen(info);                     % 添加CRC校验位到信息比特中

% ----------------------高斯近似构造极化码的信道信息-------------------------
sigma = 1 / sqrt(2 * R) * 10^(-snr / 20); % 高斯信道信噪比方差,sigma表示高斯白噪声的标准差
channels = GA(sigma, N);                  % 使用高斯近似构造信道质量
% 对 channels 按信道质量从高到低排序
[~, channel_ordered] = sort(channels, 'descend'); 
% 从 channel_ordered 中取出前K 个最可靠的信道作为信息位
info_bits = sort(channel_ordered(1 : K), 'ascend'); 

frozen_bits = ones(N, 1); 
frozen_bits(info_bits) = 0;              % 信息位为0，冻结比特为1
logic = mod(frozen_bits + 1, 2);         % 信息位为1，冻结比特为0
info_bits_logical = logical(logic);      % 将数值转化为逻辑值

% --------------------------极化码编码及调制-------------------------------
u = zeros(N, 1);
u(info_bits_logical) = info_with_crc;    % 将信息比特放入信源向量
x = polar_encoder(u);                    % 极化码编码
% BPSK调制，x = 0 被映射为 +1，x = 1 被映射为 -1。
% 0 映射为相位 0°，即 +1，1 映射为相位 180°，即 -1。
bpsk = 1 - 2 * x;                        

% ----------------------------仿真循环部分-------------------------------
blenum_sc = 0;  % SC译码误块次数初始化
blenum_scl = 0; % SCL译码误块次数初始化
blenum_cascl = 0; % CASCL译码误块次数初始化

for i_run = 1 : max_runs
    % 添加噪声并计算信道LLR
    y = awgn(bpsk, snr);                 % 添加高斯白噪声(Additive White Gaussian Noise)
    llr = 2 / sigma^2 * y;               % 计算信道LLR

    % 译码
    polar_info_sc = SC_decoder(llr, K, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec); % SC译码
    polar_info_scl = SCL_decoder(llr, L, K, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec); % SCL译码
    % polar_info_cascl = CASCL_decoder(llr, L, K, frozen_bits, det, lambda_offset, llr_layer_vec, bit_layer_vec); % CA-SCL译码

    % 统计误块次数
    if any(polar_info_sc ~= info_with_crc)
        blenum_sc = blenum_sc + 1;       % SC译码失败累计
    end
    if any(polar_info_scl ~= info_with_crc)
        blenum_scl = blenum_scl + 1;     % SCL译码失败累计
    end
    % if any(polar_info_cascl ~= info_with_crc)
    %     blenum_cascl = blenum_cascl + 1; % CA-SCL译码失败累计
    % end
end

% --------------------------统计误块率及输出结果----------------------------
bler_sc = blenum_sc / max_runs;          % SC误块率
bler_scl = blenum_scl / max_runs;        % SCL误块率
bler_cascl = blenum_cascl / max_runs;    % CA-SCL误块率

fprintf('Sim iteration running = %d\n', max_runs); % 仿真迭代运行次数
fprintf('N = %d, K = %d, L = %d\n', N, K, L);     % 输出极化码参数
fprintf('The SNR = %.1f\n', snr);                 % 输出信噪比
fprintf('The BLER of SC = %f\n', bler_sc);        % SC译码误块率
fprintf('The BLER of SCL = %f\n', bler_scl);      % SCL译码误块率
% fprintf('The BLER of CA-SCL = %f\n', bler_cascl); % CA-SCL译码误块率
