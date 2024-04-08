function para_out = para_exp_transform_simplify(para_log,is_posi_nega)
para_out = para_log;

para_out(is_posi_nega(:,1),:) = exp(para_log(is_posi_nega(:,1),:));
para_out(is_posi_nega(:,2),:) =-exp(para_log(is_posi_nega(:,2),:));
end