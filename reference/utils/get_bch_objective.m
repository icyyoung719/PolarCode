function [enc, dec] = get_bch_objective(n, t)
% 根据极化码参数n和纠错能力t生成BCH编解码器
%   n: 极化码指数（N=2^n）
%   t: BCH纠错能力
%   enc: BCH编码器对象
%   dec: BCH解码器对象

N_polar = 2^n;
m = nextpow2(N_polar + 1);  % 确保n_bch ≥ N_polar
n_bch = 2^m - 1;            % BCH码长
possible = bchnumerr(n_bch);% 获取所有可能参数组合

% 筛选满足纠错能力t的参数
valid = possible(possible(:,2) >= t, :);
if isempty(valid)
    error('No valid BCH parameters for n=%d, t=%d', n, t);
end

% 选择最大信息位长度
[~, idx] = max(valid(:,1));
k_bch = valid(idx, 1);
t_actual = valid(idx, 2);

% 创建BCH编解码器
enc = comm.BCHEncoder(n_bch, k_bch, ...
    GeneratorPolynomialSource='Auto', ...
    PuncturePatternSource='None');

dec = comm.BCHDecoder(n_bch, k_bch, ...
    GeneratorPolynomialSource='Auto', ...
    NumCorrectedErrorsOutputPort=true);
end