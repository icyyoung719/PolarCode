function result = f_l(llr_left, llr_right, u_hat_left)
    % 左分支的LLR更新规则
    result = log(exp(llr_left + llr_right) + 1) - log(exp(llr_left) + exp(llr_right));
end