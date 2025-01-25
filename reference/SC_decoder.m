function polar_info_esti = SC_decoder(llr, K, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec)
% 输入：
%   llr            - 信道 LLR 值（对信道接收值的可靠性估计）
%   K              - 信息比特的个数
%   frozen_bits    - 冻结比特位置标识，长度为 N（1 表示冻结比特，0 表示信息比特）
%   lambda_offset  - 每一层的偏移量，用于定位递归计算的范围
%   llr_layer_vec  - 每个比特对应的 LLR 层编号
%   bit_layer_vec  - 每个比特对应的译码层编号
% 输出：
%   polar_info_esti - 译码出的信息比特序列

% 参数初始化
N = length(llr);                   % 码长（必须为 2 的幂次）
n = log2(N);                       % 码的深度（n = log2(N)）
P = zeros(N - 1, 1);               % 存储内部节点的 LLR 值（不包括叶子节点）
C = zeros(N - 1, 2) - 1;           % 存储内部节点的比特值（2 列分别对应左/右路径）
polar_info_esti = zeros(K, 1);     % 存储译码得到的信息比特
cnt_K = 1;                         % 记录当前已译码出信息比特的数量

for phi = 0 : N - 1
    switch phi
    % 根据 𝜙 的位置，确定其属于哪个阶段（译码第一半部分、第二半部分或一般位置）。
        case 0 % 译码上一半的第一个
            index_1 = lambda_offset(n);  % index_1 是 4， 如果 N=8 的话
            % 这个 for 循环，计算倒数第二层的 LLR 似然比，N=8极化码，计算出来 4 个 LLR
            % P(4) <==  LLR(1) 与 LLR(5)
            % P(5) <==  LLR(2) 与 LLR(6)
            % P(6) <==  LLR(3) 与 LLR(7)
            % P(7) <==  LLR(4) 与 LLR(8)
            for beta = 0 : index_1 - 1 % use llr vector
                % L_left = sign(L_1)*sign(L_2)*min{abs(L_1),abs(L_2)}
                P(beta + index_1) =  sign(llr(beta + 1)) * sign(llr(beta + 1 + index_1)) * min(abs(llr(beta + 1)), abs(llr(beta + 1 + index_1)));
            end
            % 这个 for 循环，计算除了倒数第二层其它层的 LLR 似然比，循环两次，分别计算都输第三层(2个）和第四层（1个）的 LLR
            % i_layer = 1 时, index_1=2  index_2 = 4
            %    P(2)  <== P(4) 与 P(6)
            %    P(3)  <== P(5) 与 P(7)
            % i_layer = 0 时, index_1=1  index_2 = 2
            %    P(1)  <== P(2) 与 P(3)
            for i_layer = n - 2 : -1 : 0 % use P vector
                index_1 = lambda_offset(i_layer + 1);
                index_2 = lambda_offset(i_layer + 2);
                for beta = index_1 : index_2 - 1
                    P(beta) =  sign(P(beta + index_1)) * sign(P(beta + index_2)) * min(abs(P(beta + index_1)), abs(P(beta + index_2)));
                end
            end
        case N/2 % 译码下一半的第一个
            index_1 = lambda_offset(n);
            % 下一半的倒数第二层，肩膀上都有译码出来的比特
            for beta = 0 : index_1 - 1 % use llr vector. g function.
                % 右分支信道的 LLR：考虑已译码出的比特
                % L_right = (1 - 2*u)*L_1 + L_2, u是先前译码出的比特值
                P(beta + index_1) = (1 - 2 * C(beta + index_1, 1)) * llr(beta + 1) + llr(beta + 1 + index_1);
            end
            % 再译码向上（向左）各层的似然比，肩膀没有译码出来的比特的情况
            for i_layer = n - 2 : -1 : 0 % use P vector. f function
                % 继续向上合并 LLR
                index_1 = lambda_offset(i_layer + 1);
                index_2 = lambda_offset(i_layer + 2);
                for beta = index_1 : index_2 - 1
                    P(beta) =  sign(P(beta + index_1)) * sign(P(beta + index_2)) * min(abs(P(beta + index_1)), abs(P(beta + index_2)));
                end
            end
        otherwise
            % 从当前比特对应的 LLR 层向根节点递归计算
            llr_layer = llr_layer_vec(phi + 1);
            index_1 = lambda_offset(llr_layer + 1);
            index_2 = lambda_offset(llr_layer + 2);
            % 这个是计算前面已经有译码出来 u 的N=2极化码的 似然比
            % 使用公式 L = (1 - 2*u)*L_1 + L_2
            for beta = index_1 : index_2 - 1 % g function is first implemented.
                P(beta) = (1 - 2 * C(beta, 1)) * P(beta + index_1) + P(beta + index_2);
            end
            % 计算 N=2 极化码的没有译码出比特的情况，可能有多层
            for i_layer = llr_layer - 1 : -1 : 0 % then f function is implemented.
                index_1 = lambda_offset(i_layer + 1);
                index_2 = lambda_offset(i_layer + 2);
                % 使用公式 L = sign(L_1)*sign(L_2)*min{abs(L_1),abs(L_2)}
                for beta = index_1 : index_2 - 1
                    P(beta) =  sign(P(beta + index_1)) * sign(P(beta + index_2)) * min(abs(P(beta + index_1)), abs(P(beta + index_2)));
                end
            end
    end
    phi_mod_2 = mod(phi, 2);
% -----------------------------冻结比特与信息比特处理----------------------
    if frozen_bits(phi + 1) == 1 % 冻结比特
        C(1, 1 + phi_mod_2) = 0; % 冻结比特固定为 0
    else % 信息比特
        % 最终需要计算的值不断被归约到 P(1)，根节点
        C(1, 1 + phi_mod_2) = P(1) < 0; % LLR < 0 表示估计比特为 1
        polar_info_esti(cnt_K) = P(1) < 0; % 保存估计值
        cnt_K = cnt_K + 1;
    end
% -----------------------------更新内部节点比特值-------------------------
% 更新内部节点比特值 C 是为了完成极化码的递归回溯操作，将已译码比特的信息向上层传播。
% 极化码的译码过程类似二叉树结构，每次译码一个比特，都会影响到对应内部节点的状态，进而影响后续比特的译码。
    if phi_mod_2  == 1 && phi ~= N - 1
        bit_layer = bit_layer_vec(phi + 1); % 当前比特所在层级
        for i_layer = 0 : bit_layer - 1 % give values to the 2nd column of C
            index_1 = lambda_offset(i_layer + 1);
            index_2 = lambda_offset(i_layer + 2);
            for beta = index_1 : index_2 - 1
                % 遍历当前层的所有节点,更新左分支和右分支的比特值
                % u_left = (u_parent + u_right) mod 2
                % u_right = u_right 
                C(beta + index_1, 2) = mod(C(beta, 1) + C(beta, 2), 2);
                C(beta + index_2, 2) = C(beta, 2);
            end
        end
        index_1 = lambda_offset(bit_layer + 1);
        index_2 = lambda_offset(bit_layer + 2);
        for beta = index_1 : index_2 - 1 % give values to the 1st column of C
            C(beta + index_1, 1) = mod(C(beta, 1) + C(beta, 2), 2);
            C(beta + index_2, 1) = C(beta, 2);
        end
    end
end
end

