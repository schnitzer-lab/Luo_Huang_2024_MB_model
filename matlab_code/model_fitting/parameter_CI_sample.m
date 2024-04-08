function [x_CI, x_rand] = parameter_CI_sample(x_mu,hessian,sample_n,alpha)
% alpha is the confidence level, default value: 0.05
if ~exist('alpha','var')
    alpha = 0.05;
end
% log(likelihood) = -1/2 * Error_sq_sum
% Fisher_info = -diff(log(likelihood),2)
% Hessian = diff(Error_sq_sum, 2)
Fisher_info = hessian/2;
x_sigma_squre = Fisher_info^-1;

x_CI = x_mu - norminv(alpha/2)*sqrt(diag(x_sigma_squre))*[-1 1];
x_rand = mvnrnd(x_mu',x_sigma_squre,sample_n);
end
