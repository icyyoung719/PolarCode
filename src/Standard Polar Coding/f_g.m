function result = f_g(llr_left, llr_right, u_hat_left)
    % 右分支的LLR更新规则
    result = llr_right + (-1).^u_hat_left .* llr_left;
end
