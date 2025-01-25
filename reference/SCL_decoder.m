function polar_info_esti = SCL_decoder(llr, L, K, frozen_bits, lambda_offset, llr_layer_vec, bit_layer_vec)
% 输入：
%   llr            - 信道 LLR 值
%   L              - 候选路径的最大数量
%   K              - 信息比特的个数
%   frozen_bits    - 冻结比特位置标识，长度为 N（1 表示冻结比特，0 表示信息比特）
%   lambda_offset  - 每一层的偏移量，用于定位递归计算的范围
%   llr_layer_vec  - 每个比特对应的 LLR 层编号
%   bit_layer_vec  - 每个比特对应的译码层编号
% 输出：
%   polar_info_esti - 译码出的信息比特序列

% -----------------------------初始化--------------------------------------
N = length(llr);                  % 码长（必须为 2 的幂次）
m = log2(N);                      % 码的深度（m = log2(N)）

% 内存分配
lazy_copy = zeros(m, L);          % 懒拷贝数组，用于记录数据来源路径，避免数据复制
P = zeros(N - 1, L);              % 存储内部节点的 LLR 值（不包括叶子节点）
C = zeros(N - 1, 2 * L) - 1;      % 存储内部节点的比特值（2 列分别对应左/右路径）
u = zeros(K, L);                  % 存储已经译出来的信息比特
PM = zeros(L, 1);                 % 路径度量
activepath = zeros(L, 1);         % 活跃路径指示器（1 表示活跃，0 表示非活跃）
polar_info_esti = zeros(K, 1);    % 译码得到的信息比特序列
cnt_u = 1;                        % 信息比特计数器

activepath(1) = 1;                % 初始化路径 0 为活跃路径
lazy_copy(:, 1) = 1;              % 懒拷贝初始化

% -----------------------------开始译码--------------------------------------
%default: in the case of path clone, the origianl path always corresponds to bit 0, while the new path bit 1.
for phi = 0 : N - 1
    layer = llr_layer_vec(phi + 1); % 当前比特所在的 LLR 层
    phi_mod_2 = mod(phi, 2);        % 当前比特位置奇偶性（用于部分求和更新）

    % 更新每个活跃路径上的 LLR
    for l_index = 1 : L
        if activepath(l_index) == 0
            continue;
        end
        switch phi
            % 根据 𝜙 的位置，确定其属于哪个阶段
            case 0
                index_1 = lambda_offset(m);
                for beta = 0 : index_1 - 1 % 计算倒数第二层的 LLR 似然比
                    P(beta + index_1, l_index) = sign(llr(beta + 1)) * sign(llr(beta + index_1 + 1)) * min(abs(llr(beta + 1)), abs(llr(beta + index_1 + 1)));
                end
                
                for i_layer = m - 2 : -1 : 0 % 计算除了倒数第二层其它层的 LLR 似然比
                    index_1 = lambda_offset(i_layer + 1);
                    index_2 = lambda_offset(i_layer + 2);
                    for beta = 0 : index_1 - 1
                        P(beta + index_1, l_index) = sign(P(beta + index_2, l_index)) *...
                            sign(P(beta + index_1 + index_2, l_index)) * min(abs(P(beta + index_2, l_index)), abs(P(beta + index_1 + index_2, l_index)));
                    end
                end
            case N/2
                index_1 = lambda_offset(m);
                for beta = 0 : index_1 - 1
                    x_tmp = C(beta + index_1, 2 * l_index - 1);
                    P(beta + index_1, l_index) = (1 - 2 * x_tmp) * llr(beta + 1) + llr(beta + 1 + index_1);
                end
                for i_layer = m - 2 : -1 : 0
                    index_1 = lambda_offset(i_layer + 1);
                    index_2 = lambda_offset(i_layer + 2);
                    for beta = 0 : index_1 - 1
                        P(beta + index_1, l_index) = sign(P(beta + index_2, l_index)) *...
                            sign(P(beta + index_1 + index_2, l_index)) * min(abs(P(beta + index_2, l_index)), abs(P(beta + index_1 + index_2, l_index)));
                    end
                end
            otherwise
                index_1 = lambda_offset(layer + 1);
                index_2 = lambda_offset(layer + 2);
                for beta = 0 : index_1 - 1
                    P(beta + index_1, l_index) = (1 - 2 * C(beta + index_1, 2 * l_index - 1)) * P(beta + index_2, lazy_copy(layer + 2, l_index)) +...
                        P(beta + index_1 + index_2, lazy_copy(layer + 2, l_index));
                end
                for i_layer = layer - 1 : -1 : 0
                    index_1 = lambda_offset(i_layer + 1);
                    index_2 = lambda_offset(i_layer + 2);
                    for beta = 0 : index_1 - 1
                        P(beta + index_1, l_index) = sign(P(beta + index_2, l_index)) *...
                            sign(P(beta + index_1 + index_2, l_index)) * min(abs(P(beta + index_2, l_index)),...
                            abs(P(beta + index_1 + index_2, l_index)));
                    end
                end
        end
    end
% -----------------------------冻结比特与信息比特处理----------------------
    if frozen_bits(phi + 1) == 0 % 信息比特
        PM_pair = realmax * ones(2, L); % 初始化路径度量，每一列表示每条路径
        % 更新每个路径上的度量值
        for l_index = 1 : L
            if activepath(l_index) == 0
                continue;
            end
            if P(1, l_index) >= 0 % +绝对值
                PM_pair(1, l_index) = PM(l_index);
                PM_pair(2, l_index) = PM(l_index) + P(1, l_index);
            else
                PM_pair(1, l_index) = PM(l_index) - P(1, l_index);
                PM_pair(2, l_index) = PM(l_index);
            end
        end
        % 排序路径度量，选择最优路径
        middle = min(2 * sum(activepath), L);
        PM_sort = sort(PM_pair(:));
        PM_cv = PM_sort(middle);
        % compare 中为1的留下，为0的丢弃
        compare = PM_pair <= PM_cv; 

        % 更新路径状态
        kill_index = zeros(L, 1); % 用于记录被移除路径的索引
        kill_cnt = 0; % 被移除路径计数器
        for i = 1 : L
            if (compare(1, i) == 0)&&(compare(2, i) == 0)%which indicates that this path should be killed
                activepath(i) = 0;
                kill_cnt = kill_cnt + 1;%push stack
                kill_index(kill_cnt) = i;
            end
        end
% -----------------------------路径克隆与选择------------------------------
        for l_index = 1 : L
            if activepath(l_index) == 0
                continue;
            end
            path_state = compare(1, l_index) * 2 + compare(2, l_index);
            % path_state: 当前路径 l_index 下新长出的两个路径（两种可能，当前译码的比特取0或者1）
            %   1：保留译码为1的路径
            %   2：保留译码为0的路径
            %   3：保留两条路径
            %   0：不保留
            switch path_state%path_state can equal to 0, but in this case we do no operation.
                case 1
                    u(cnt_u, l_index) = 1;
                    C(1, 2 * l_index - 1 + phi_mod_2) = 1;
                    PM(l_index) = PM_pair(2, l_index);
                case 2
                    u(cnt_u, l_index) = 0;
                    C(1, 2 * l_index - 1 + phi_mod_2) = 0;
                    PM(l_index) = PM_pair(1, l_index);
                case 3
                    % 克隆路径
                    index = kill_index(kill_cnt);  % 从最后一个 kill 的路径
                    kill_cnt = kill_cnt - 1; % pop stack
                    activepath(index) = 1;
                    % 懒拷贝并更新比特
                    lazy_copy(:, index) = lazy_copy(:, l_index);
                    u(:, index) = u(:, l_index);
                    u(cnt_u, l_index) = 0;
                    u(cnt_u, index) = 1;
                    C(1, 2 * l_index - 1 + phi_mod_2) = 0;
                    C(1, 2 * index - 1 + phi_mod_2) = 1;
                    PM(l_index) = PM_pair(1, l_index);
                    PM(index) = PM_pair(2, l_index);
            end
        end
        cnt_u = cnt_u + 1;
    else % 冻结比特
        for l_index = 1 : L
            if activepath(l_index) == 0
                continue;
            end
            if P(1, l_index) < 0
                PM(l_index) = PM(l_index) - P(1, l_index);
            end
            if phi_mod_2 == 0
                C(1, 2 * l_index - 1) = 0;
            else
                C(1, 2 * l_index) = 0;
            end 
        end
    end 
    
% -----------------------------更新内部节点比特值-------------------------
    for l_index = 1 : L%partial-sum return
        if activepath(l_index) == 0
            continue
        end
        if (phi_mod_2  == 1) && (phi ~= N - 1)
            layer = bit_layer_vec(phi + 1);
            for i_layer = 0 : layer - 1
                index_1 = lambda_offset(i_layer + 1);
                index_2 = lambda_offset(i_layer + 2);
                for beta = index_1 : 2 * index_1 - 1
                    C(beta + index_1, 2 * l_index) = mod(C(beta, 2 *  lazy_copy(i_layer + 1, l_index) - 1) + C(beta, 2 * l_index), 2);%Left Column lazy copy
                    C(beta + index_2, 2 * l_index) = C(beta, 2 * l_index);   
                end
            end
            index_1 = lambda_offset(layer + 1);
            index_2 = lambda_offset(layer + 2); 
            for beta = index_1 : 2 * index_1 - 1
                C(beta + index_1, 2 * l_index - 1) = mod(C(beta, 2 * lazy_copy(layer + 1, l_index) - 1) + C(beta, 2 * l_index), 2);%Left Column lazy copy
                C(beta + index_2, 2 * l_index - 1) = C(beta, 2 * l_index);
            end 
        end
    end
    % lazy copy
    if phi < N - 1
        for i_layer = 1 : llr_layer_vec(phi + 2) + 1
            for l_index = 1 : L
                lazy_copy(i_layer, l_index) = l_index;
            end
        end
    end
end
% ---------------------------------输出------------------------------------
% 路径选择：选择路径度量最小的路径
[~, min_index] = min(PM); %输出按升序排列的PM值的索引
polar_info_esti=u(:,min_index);
end
