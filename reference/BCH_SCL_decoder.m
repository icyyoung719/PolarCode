function polar_info_esti = BCH_SCL_decoder(llr, L, K, frozen_bits, ...
   dec_bch,lambda_offset, llr_layer_vec, bit_layer_vec)
% 极化码CASCL译码器（LLR-based，单函数实现）
% 参数：信道LLR、列表长度L、信息位数K、冻结位标识、CRC校验函数、层偏移量等

N = length(llr);          % 码长
m = log2(N);              % 码层数
lazy_copy = zeros(m, L);  % 惰性拷贝记录
P = zeros(N-1, L);        % 节点LLR存储
C = zeros(N-1, 2*L);      % 节点比特存储
u = zeros(K, L);          % 信息位存储（含CRC）
PM = zeros(L, 1);         % 路径度量
activepath = zeros(L,1);  % 活跃路径标记
cnt_u = 1;                % 信息位计数器

activepath(1) = 1;        % 初始化第一条路径
lazy_copy(:,1) = 1;

% 主解码循环
for phi = 0:N-1
    layer = llr_layer_vec(phi+1);
    phi_mod_2 = mod(phi,2);
    
    for l_index = 1:L
        if ~activepath(l_index), continue; end
        
        % 特殊节点处理（u_0和u_N/2）
        switch phi
            case 0
                index_1 = lambda_offset(m);
                for beta = 0:index_1-1
                    P(beta+index_1, l_index) = sign(llr(beta+1)) * sign(llr(beta+index_1+1)) * ...
                        min(abs(llr(beta+1)), abs(llr(beta+index_1+1)));
                end
                for i_layer = m-2:-1:0
                    index_1 = lambda_offset(i_layer+1);
                    index_2 = lambda_offset(i_layer+2);
                    for beta = 0:index_1-1
                        P(beta+index_1, l_index) = sign(P(beta+index_2, l_index)) * ...
                            sign(P(beta+index_1+index_2, l_index)) * ...
                            min(abs(P(beta+index_2, l_index)), abs(P(beta+index_1+index_2, l_index)));
                    end
                end
            case N/2
                index_1 = lambda_offset(m);
                for beta = 0:index_1-1
                    x_tmp = C(beta+index_1, 2*l_index-1);
                    P(beta+index_1, l_index) = (1-2*x_tmp)*llr(beta+1) + llr(beta+1+index_1);
                end
                for i_layer = m-2:-1:0
                    index_1 = lambda_offset(i_layer+1);
                    index_2 = lambda_offset(i_layer+2);
                    for beta = 0:index_1-1
                        P(beta+index_1, l_index) = sign(P(beta+index_2, l_index)) * ...
                            sign(P(beta+index_1+index_2, l_index)) * ...
                            min(abs(P(beta+index_2, l_index)), abs(P(beta+index_1+index_2, l_index)));
                    end
                end
            otherwise
                index_1 = lambda_offset(layer+1);
                index_2 = lambda_offset(layer+2);
                for beta = 0:index_1-1
                    P(beta+index_1, l_index) = (1-2*C(beta+index_1, 2*l_index-1)) * ...
                        P(beta+index_2, lazy_copy(layer+2, l_index)) + ...
                        P(beta+index_1+index_2, lazy_copy(layer+2, l_index));
                end
                for i_layer = layer-1:-1:0
                    index_1 = lambda_offset(i_layer+1);
                    index_2 = lambda_offset(i_layer+2);
                    for beta = 0:index_1-1
                        P(beta+index_1, l_index) = sign(P(beta+index_2, l_index)) * ...
                            sign(P(beta+index_1+index_2, l_index)) * ...
                            min(abs(P(beta+index_2, l_index)), abs(P(beta+index_1+index_2, l_index)));
                    end
                end
        end
    end
    
    % 冻结位/信息位处理
    if frozen_bits(phi+1) == 0  % 信息位
        PM_pair = realmax*ones(2,L);
        for l_index = 1:L
            if ~activepath(l_index), continue; end
            if P(1,l_index) >= 0
                PM_pair(:,l_index) = [PM(l_index); PM(l_index)+P(1,l_index)];
            else
                PM_pair(:,l_index) = [PM(l_index)-P(1,l_index); PM(l_index)];
            end
        end
        
        % 路径筛选与克隆
        middle = min(2*sum(activepath), L);
        PM_sort = sort(PM_pair(:));
        PM_cv = PM_sort(middle);
        compare = PM_pair <= PM_cv;
        
        kill_index = zeros(L,1);
        kill_cnt = 0;
        for i = 1:L
            if all(compare(:,i) == 0)
                activepath(i) = 0;
                kill_cnt = kill_cnt + 1;
                kill_index(kill_cnt) = i;
            end
        end
        
        for l_index = 1:L
            if ~activepath(l_index), continue; end
            path_state = compare(1,l_index)*2 + compare(2,l_index);
            
            switch path_state
                case 1  % 保留bit=1路径
                    u(cnt_u, l_index) = 1;
                    C(1, 2*l_index-1+phi_mod_2) = 1;
                    PM(l_index) = PM_pair(2,l_index);
                case 2  % 保留bit=0路径
                    u(cnt_u, l_index) = 0;
                    C(1, 2*l_index-1+phi_mod_2) = 0;
                    PM(l_index) = PM_pair(1,l_index);
                case 3  % 克隆新路径
                    index = kill_index(kill_cnt);
                    kill_cnt = kill_cnt - 1;
                    activepath(index) = 1;
                    
                    % 惰性拷贝
                    lazy_copy(:,index) = lazy_copy(:,l_index);
                    u(:,index) = u(:,l_index);
                    u(cnt_u, l_index) = 0;
                    u(cnt_u, index) = 1;
                    C(1, 2*l_index-1+phi_mod_2) = 0;
                    C(1, 2*index-1+phi_mod_2) = 1;
                    PM(l_index) = PM_pair(1,l_index);
                    PM(index) = PM_pair(2,l_index);
            end
        end
        cnt_u = cnt_u + 1;
    else  % 冻结位处理
        for l_index = 1:L
            if ~activepath(l_index), continue; end
            if P(1,l_index) < 0
                PM(l_index) = PM(l_index) - P(1,l_index);
            end
            if phi_mod_2 == 0
                C(1, 2*l_index-1) = 0;
            else
                C(1, 2*l_index) = 0;
            end 
        end
    end
    
    % 部分和更新
    for l_index = 1:L
        if ~activepath(l_index) || phi == N-1 || ~phi_mod_2, continue; end
        layer = bit_layer_vec(phi+1);
        
        for i_layer = 0:layer-1
            index_1 = lambda_offset(i_layer+1);
            index_2 = lambda_offset(i_layer+2);
            for beta = index_1:2*index_1-1
                C(beta+index_1, 2*l_index) = mod(C(beta, 2*lazy_copy(i_layer+1,l_index)-1) + ...
                    C(beta, 2*l_index), 2);
                C(beta+index_2, 2*l_index) = C(beta, 2*l_index);
            end
        end
        
        index_1 = lambda_offset(layer+1);
        index_2 = lambda_offset(layer+2);
        for beta = index_1:2*index_1-1
            C(beta+index_1, 2*l_index-1) = mod(C(beta, 2*lazy_copy(layer+1,l_index)-1) + ...
                C(beta, 2*l_index), 2);
            C(beta+index_2, 2*l_index-1) = C(beta, 2*l_index);
        end 
    end
    
    % 惰性拷贝重置
    if phi < N-1
        for i_layer = 1:llr_layer_vec(phi+2)+1
            lazy_copy(i_layer,:) = 1:L;
        end
    end
end

% CRC校验路径选择
[~, path_ordered] = sort(PM);
for l_index = 1:L
    path_num = path_ordered(l_index);
    info_with_bch = u(:, path_num);
    [~,err] = dec_bch(info_with_bch);
    if err == 0
        polar_info_esti = u(:, path_num);
        break;
    elseif l_index == L
        polar_info_esti = u(:, path_ordered(1));
    end
end
end