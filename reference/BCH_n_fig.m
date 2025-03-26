% 参数配置
max_runs = 100;  % 仿真次数
L = 16;           % SCL列表长度
t_bch = 2;        % BCH纠错能力

% =============================================================
% 图1：不同SNR下的BLER性能（固定n=8）


snr = 2.0;              % 信噪比（dB）
n_range = 6:10;  % n从4到10（N=16到1024）
bler_sc_n = zeros(size(n_range));
bler_scl_n = zeros(size(n_range));
bler_bchscl_n = zeros(size(n_range));

for idx = 1:length(n_range)
    n = n_range(idx);
    N = 2^n;
    msgbit_len = 26 * (2^(n-6)); % 信息位长度
    R = msgbit_len / N;  % 码率
    K = msgbit_len;      % 编码后的总长度位

    fprintf('Processing N=%d dB\n', n);
    
    % 初始化误块计数器
    blenum_sc = 0;
    blenum_scl = 0;
    blenum_bchscl = 0;

    % 动态生成BCH参数
    S = 6;         % Shortened message length
    gp = 'x5+x2+1'; % Generator polynomial
    enc_bch = comm.BCHEncoder(31, 26, gp);
    dec_bch = comm.BCHDecoder(31, 26, gp);
    n_bch = enc_bch.CodewordLength;
    k_bch = enc_bch.MessageLength;

    bch_total_len = ceil(msgbit_len / k_bch) * n_bch; % 总码字长度
    K = bch_total_len;                                % 编码后的总长度位
    bch_parity_len = n_bch - k_bch;                  % BCH校验位长度

    % 预处理（仅需执行一次）
    lambda_offset = 2.^(0:log2(N));  % 分段向量，表示每一层的信道分段长度
    llr_layer_vec = get_llr_layer(N);  % LLR计算的实际执行层数向量
    bit_layer_vec = get_bit_layer(N);  % 比特值返回时实际执行层数向量

    % 信道构造
    sigma = 1 / sqrt(2 * R) * 10^(-snr / 20);
    channels = GA(sigma, N);
    [~, channel_ordered] = sort(channels, 'descend');
    info_bits = sort(channel_ordered(1:K), 'ascend');
    frozen_bits = ones(N, 1);
    frozen_bits(info_bits) = 0;
    logic = mod(frozen_bits + 1, 2);
    info_bits_logical = logical(logic);

    for run = 1:max_runs
        % 生成信息比特并添加BCH校验
        info = randi([0 1], msgbit_len, 1); % 生成k_bch位信息
        info_with_bch = enc_bch(info);     % BCH编码
        info_with_bch = info_with_bch(:);  % 转为列向量

        % 极化码编码与调制
        u = zeros(N, 1);
        u(info_bits_logical) = info_with_bch;
        x = polar_encoder(u);
        bpsk = 1 - 2 * x;

        % 加入高斯噪声
        y = awgn(bpsk, snr);
        llr = 2 / sigma^2 * y;

        % 译码
        polar_info_sc = SC_decoder(llr, K, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec);
        polar_info_scl = SCL_decoder(llr, L, K, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec);
        polar_info_bchscl = BCH_SCL_decoder(llr, L, K, frozen_bits, dec_bch, lambda_offset, llr_layer_vec, bit_layer_vec);

        % 统计误块
        if any(polar_info_sc ~= info_with_bch)
            blenum_sc = blenum_sc + 1;
        end
        if any(polar_info_scl ~= info_with_bch)
            blenum_scl = blenum_scl + 1;
        end
        if any(polar_info_bchscl ~= info_with_bch)
            blenum_bchscl = blenum_bchscl + 1;
        end
    end

    % 计算BLER
    bler_sc_n(idx) = blenum_sc / max_runs;
    bler_scl_n(idx) = blenum_scl / max_runs;
    bler_bchscl_n(idx) = blenum_bchscl / max_runs;
end

% 绘制n-BLER曲线
figure;
semilogy(n_range, bler_sc_n, 'bo-', 'LineWidth', 1.5);
hold on;
semilogy(n_range, bler_scl_n, 'rs-', 'LineWidth', 1.5);
semilogy(n_range, bler_bchscl_n, 'gv-', 'LineWidth', 1.5);
xlabel('n (N=2^n)');
ylabel('BLER');
title('BLER Performance vs Code Length (SNR=2dB)');
legend('SC', 'SCL (L=16)', 'BCH-SCL (L=16)');
grid on;
axis([6 10 1e-4 1]);