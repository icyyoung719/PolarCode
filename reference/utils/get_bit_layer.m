function layer_vec = get_bit_layer(N)
    % 初始化一个长度为 N 的向量，用于存储每个比特的回传层级
    layer_vec = zeros(N, 1);  
    
    % 遍历所有比特索引，从 0 到 N-1
    for phi = 0 : N - 1  
        % 当前比特索引整除 2
        psi = floor(phi / 2);  
        % 初始化层数为 0
        layer = 0;  
        
        % 当 psi 的最低有效位（LSB）为 1 时，继续计算层数
        while(mod(psi, 2) == 1)  
            psi = floor(psi / 2);  % 整除 2
            layer = layer + 1;    % 层数加 1
        end
        
        % 将计算得到的层数存入对应索引位置
        layer_vec(phi + 1) = layer;  
    end
end
