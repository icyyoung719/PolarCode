function x = phi_inverse(y)
%部分用闭合表达式，部分用数值解法，速度进一步提升！
if (y <= 1.0221) && (y >= 0.0388)
    x = ((0.0218 - log(y))/0.4527)^(1/0.86);
else
    x0 = 0.0388;
    x1 = x0 - (phi(x0) - y)/derivative_phi(x0);
    delta = abs(x1 - x0);
    epsilon = 1e-3;
    
    while(delta >= epsilon)
        x0 = x1;
        x1 = x1 - (phi(x1) - y)/derivative_phi(x1);
        %当x1过大，放宽epsilon
        if x1 > 1e2
            epsilon = 10;
        end       
        delta = abs(x1 - x0);
    end
    x = x1;
end
end
