function [fig_h,LL_mat] = plot_likelihood_ratio(x_mu,x_CI,fun_LL_ratio,alpha)
% alpha is the confidence level, default value: 0.05
if ~exist('alpha','var')
    alpha = 0.05;
end

fig_h = figure;
x_n = size(x_mu,1);
axis_h = zeros(x_n);
x_change_n = 31;
LL_mat = zeros(x_n,x_change_n);
for x_i = 1:x_n
    axis_h(x_i) = subplot(ceil(sqrt(x_n)),ceil(sqrt(x_n)),x_i); hold on
    x_mu_curr = x_mu;
    
    x_change_list = linspace(x_CI(x_i,1),x_CI(x_i,2),x_change_n);
    for x_change_i = 1:x_change_n
        x_mu_curr(x_i) = x_change_list(x_change_i);
        LL_mat(x_i,x_change_i) = fun_LL_ratio(x_mu_curr);
    end
    
    plot(x_change_list,exp(LL_mat(x_i,:)))
    y = normpdf(linspace(-norminv(alpha/2),norminv(alpha/2),x_change_n))/normpdf(0);
    plot(x_change_list,y)
    
    set(gca,'XLim',x_CI(x_i,:))
    
end
end
%{

%%

ci = 2*sqrt(diag(x_sigma_squre));

%log_pdf_posterior = @(x)-0.5*Error_sq_sum(x') ...
%                    +log(double(all(x'>para_lu(:,1)) ...
%                              & all(x'<para_lu(:,2))));
                

log_prop_pdf = @(x,y)log(mvnpdf(x,x_mu',Fisher_info^-1));
prop_rnd = @(x)mvnrnd(x_mu',Fisher_info^-1);
%prod_rnd = @(x)mvnrnd_in_range(x_lu',Fisher_info^-1,para_lu,1);

%%
%factor_list = -1:0.1:2;
%factor_list =0.5:0.05:1.5;
%Fisher_info = Hessian_accu/2;
para_std = sqrt(diag(Fisher_info^-1));
%para_std(14) = para_std(14)*10
%para_std([16 17 19]) = para_std([16 17 19])*10;
factor_list =(-2:0.1:2);

logp_list = zeros(length(factor_list),2);

fig_h = figure;
axis_h = zeros(size(x_mu));
for theta_i = 1:length(x_mu)
    axis_h(theta_i) = subplot(4,6,theta_i); hold on
    %x_lu_curr = x_lu;
    x_lu_curr = x_mu;
    for f_i = 1:length(factor_list)
        %x_lu_curr(theta_i) = x_mu(theta_i)*factor_list(f_i);
        x_lu_curr(theta_i) = x_mu(theta_i)+ para_std(theta_i)*factor_list(f_i);

        %logp_list(f_i) = log_pdf_posterior(x_lu_curr');
        logp_list(f_i,1) = log_pdf_posterior_expx(x_lu_curr');
        logp_list(f_i,2) = log_prop_pdf(x_lu_curr',0);
    end
    logp_list(:,2) = logp_list(:,2)-logp_list((length(factor_list)+1)/2,2);
    plot(factor_list,exp(logp_list))
    %plot(factor_list,logp_list)
end
%%
%factor_list = -1:0.1:2;
%factor_list =0.5:0.05:1.5;

%Fisher_info = Hessian_accu/2;
%para_std = sqrt(diag(Fisher_info^-1));
%para_std([16 17 19]) = para_std([16 17 19])*10;
para_cov = Fisher_info^-1;
%para_cov([16 17 19],:) = para_cov([16 17 19],:)*7;
%para_cov(:,[16 17 19]) = para_cov(:,[16 17 19])*7;
[V_cov,D_cov] = eig(para_cov);
%[V_cov,D_cov] = eig(Fisher_info);
factor_list =(-2:0.1:2);

logp_list = zeros(length(factor_list),2);

figure, 
for theta_i = 1:length(x_mu)
    subplot(4,6,theta_i)
    %x_lu_curr = x_lu;
    x_lu_curr = x_mu;
    for f_i = 1:length(factor_list)
        %x_lu_curr(theta_i) = x_mu(theta_i)*factor_list(f_i);
        x_lu_curr = x_mu + V_cov(:,theta_i)*(factor_list(f_i)*sqrt(D_cov(theta_i,theta_i)));

        %logp_list(f_i) = log_pdf_posterior(x_lu_curr');
        logp_list(f_i,1) = log_pdf_posterior_expx(x_lu_curr');
        logp_list(f_i,2) = log_prop_pdf(x_lu_curr',0);
    end
    logp_list(:,2) = logp_list(:,2)-logp_list((length(factor_list)+1)/2,2);
    plot(factor_list,exp(logp_list))
    %plot(factor_list,logp_list)
    title(D_cov(theta_i,theta_i))
end
%%
%para_rand = mvnrnd_in_range(x_lu',Fisher_info^-1,para_lu,1000);
%{
para_rand = mvnrnd(x_mu',Fisher_info^-1,300);
logp_rand_list = zeros(size(para_rand,1),2);
factor_rand = zeros(size(para_rand,1),1);
for f_i = 1:length(logp_rand_list)
    x_lu_curr = para_rand(f_i,:);
    factor_rand(f_i) = (x_lu_curr-x_mu')*Fisher_info*(x_lu_curr-x_mu')';
    logp_rand_list(f_i,2) = log_prop_pdf(x_lu_curr);
    logp_rand_list(f_i,1) = log_pdf_posterior_expx(x_lu_curr);
end
figure,plot(factor_rand,logp_rand_list,'.')
set(gca,'YScale','Log')
%}


end
%}