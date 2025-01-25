function x = phi_inverse(y)
% 计算给定y值对应的逆函数值
%       y - phi函数的输出值（需要满足 y <= 1.0221 且 y >= 0.0388）
%       x - phi函数的逆值（即求解对应y值的x值）

% 如果y的值处于闭合表达式适用范围内，直接使用闭合表达式进行计算
if (y <= 1.0221) && (y >= 0.0388)
    x = ((0.0218 - log(y))/0.4527)^(1/0.86);
else
    % 否则，使用数值方法进行求解
    
    % 初始猜测值
    x0 = 0.0388;  
    % 使用牛顿法迭代求解
    % x1 是基于初始猜测x0的迭代解
    x1 = x0 - (phi(x0) - y) / derivative_phi(x0);
    
    delta = abs(x1 - x0);
    % 设置精度阈值
    epsilon = 1e-3;
    
    % 使用牛顿法迭代直到满足精度要求
    while delta >= epsilon
        x0 = x1; 
        x1 = x1 - (phi(x1) - y) / derivative_phi(x1);
                % 当x1值过大时，放宽精度要求
        if x1 > 1e2
            epsilon = 10;  % 放宽精度限制
        end
        
        % 更新误差
        delta = abs(x1 - x0);
    end
    
    x = x1;
end
end
