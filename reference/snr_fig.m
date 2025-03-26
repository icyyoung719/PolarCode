% 参数配置
max_runs = 1000;  % 仿真次数
L = 16;           % SCL列表长度
crc_len = 4;      % CRC校验位长度

% =============================================================
% 图1：不同SNR下的BLER性能（固定n=8）
n = 8;
N = 2^n;
R = 0.40625;
msgbit_len = N*R;
K = msgbit_len + crc_len;

snr_range = 0:0.3:3;  % SNR范围0-3dB，步长0.3dB
bler_sc = zeros(size(snr_range));
bler_scl = zeros(size(snr_range));
bler_cascl = zeros(size(snr_range));

for i = 1:length(snr_range)
    snr = snr_range(i);
    fprintf('Processing SNR=%.1f dB\n', snr);
    
    % 初始化误块计数器
    blenum_sc = 0;
    blenum_scl = 0;
    blenum_cascl = 0;
    
    % 预处理（仅需执行一次）
    [gen, det] = get_crc_objective(crc_len);
    lambda_offset = 2.^(0:log2(N));
    llr_layer_vec = get_llr_layer(N);
    bit_layer_vec = get_bit_layer(N);
    
    % 信道构造
    sigma = 1/sqrt(2*R)*10^(-snr/20);
    channels = GA(sigma, N);
    [~, channel_ordered] = sort(channels, 'descend');
    info_bits = sort(channel_ordered(1:K), 'ascend');
    frozen_bits = ones(N,1);
    frozen_bits(info_bits) = 0;
    logic = mod(frozen_bits+1,2);
    info_bits_logical = logical(logic);
    
    for run = 1:max_runs
        % 生成信息比特并添加CRC
        info = randsrc(msgbit_len,1,[0 1;0.5 0.5]);
        info_with_crc = gen(info);
        
        % 极化码编码与调制
        u = zeros(N,1);
        u(info_bits_logical) = info_with_crc;
        x = polar_encoder(u);
        bpsk = 1 - 2*x;
        
        % 加入高斯噪声
        y = awgn(bpsk, snr);
        llr = 2/sigma^2 * y;
        
        % 译码
        polar_info_sc = SC_decoder(llr, K, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec);
        polar_info_scl = SCL_decoder(llr, L, K, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec);
        polar_info_cascl = CASCL_decoder(llr, L, K, frozen_bits, det, lambda_offset, llr_layer_vec, bit_layer_vec);
        
        % 统计误块
        if any(polar_info_sc ~= info_with_crc)
            blenum_sc = blenum_sc + 1;
        end
        if any(polar_info_scl ~= info_with_crc)
            blenum_scl = blenum_scl + 1;
        end
        if any(polar_info_cascl ~= info_with_crc)
            blenum_cascl = blenum_cascl + 1;
        end
    end
    
    % 计算BLER
    bler_sc(i) = blenum_sc/max_runs;
    bler_scl(i) = blenum_scl/max_runs;
    bler_cascl(i) = blenum_cascl/max_runs;
end

% 绘制SNR-BLER曲线
figure;
semilogy(snr_range, bler_sc, 'bo-', 'LineWidth', 1.5);
hold on;
semilogy(snr_range, bler_scl, 'rs-', 'LineWidth', 1.5);
semilogy(snr_range, bler_cascl, 'gv-', 'LineWidth', 1.5);
xlabel('SNR (dB)');
ylabel('BLER');
title('BLER Performance vs SNR (n=8)');
legend('SC', 'SCL (L=16)', 'CA-SCL (L=16)');
grid on;
axis([0 5 1e-4 1]);