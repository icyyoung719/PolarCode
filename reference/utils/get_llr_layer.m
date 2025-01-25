function layer_vec = get_llr_layer(N)
    % 初始化一个长度为 N 的向量，用于存储每个比特的 LLR 层级
    layer_vec = zeros(N, 1);  
    
    % 遍历所有比特索引，从 1 到 N-1
    for phi = 1 : N - 1  
        % 当前比特索引值
        psi = phi;  
        % 初始化层数为 0
        layer = 0;  
        
        % 当 psi 能够被 2 整除时，说明可以进一步计算层级
        while(mod(psi, 2) == 0)  
            psi = floor(psi / 2);  % 整除 2
            layer = layer + 1;    % 层数加 1
        end
        
        % 将计算得到的层数存入对应索引位置
        layer_vec(phi + 1) = layer;  
    end
end
