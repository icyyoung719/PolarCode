function x = polar_encoder(u)
%   u - 输入的比特向量（长度为 N，包含信息位和冻结位）
%   x - 编码后的极化码比特向量（长度为 N）


%   编码通过矩阵乘法实现，编码矩阵由 Kronecker 乘积生成。

N = length(u); 

GN = get_GN(N); % GN为生成矩阵

Y = u' * GN; % 进行矩阵乘法，计算编码后的比特向量

% 对结果进行 mod 2 运算，确保结果是二进制比特向量（每位取模 2）
x = mod(Y', 2); 

end
