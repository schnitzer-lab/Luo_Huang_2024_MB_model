function para_out = para_exp_transform_simplify(para_log,is_posi_nega)
%{
Because some parameters are set to be positive or negative,
we transformed them into log scale. 
for example, if p_1>0, we scale it into log(p_1); 
             if p_1<0, we scale it into log(-p_1).
This function transform the log scaled parameters back into linear scale.
%}

para_out = para_log;

para_out(is_posi_nega(:,1),:) = exp(para_log(is_posi_nega(:,1),:));
para_out(is_posi_nega(:,2),:) =-exp(para_log(is_posi_nega(:,2),:));
end