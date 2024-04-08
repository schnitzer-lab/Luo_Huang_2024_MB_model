%clc,
clear
addpath('../../steady_state_model_fitting/nonlinear')
addpath('../../steady_state_model_fitting')
addpath('..')
%%
for is_simplified_model = [true false]
is_keep_0_weight = true;false;
is_use_ga = false;
%% load original data
%data_ori_filename = fullfile('..','..','data_Cheng','Cheng_24hr_data');
data_ori_filename = fullfile('..','..','data_original','Cheng_24hr_data (revision)');

data_ori = cat(4, xlsread(data_ori_filename,'ACVvsETA'), ...
                  xlsread(data_ori_filename,'OCTvsBEN'));

%{
Dimension 1: PPL1-v1pedc, PPL1-a'2a2, PPL1-a3, MBON-v1pedc, MBON-a2sc, MBON-a3
Dimension 2: CS+, CS-
Dimension 3: Pre-traing, after 3x training, after 6x training, 1hr memory, 3hr, 24hr
Dimension 4: ACV vs ETA, OCT vs BEN
%}
Dx_mean = permute(cat(3,data_ori(1:3:end,1:6,:,:),data_ori(1:3:end,9:14,:,:)),[1 3 2 4]);
Dx_SEM  = permute(cat(3,data_ori(2:3:end,1:6,:,:),data_ori(2:3:end,9:14,:,:)),[1 3 2 4]);

if is_simplified_model
    % not include PPL1-a'2a2 MBON-a2sc data
    Dx_mean([2 5],:,:,:) = nan;
    Dx_SEM([2 5],:,:,:) = nan;
end
%{
% not include the 24hr data
Dx_mean = Dx_mean(:,:,1:end-1,:);
Dx_SEM = Dx_SEM(:,:,1:end-1,:);
%}
%% set parameter range

n_odor = 2;
n_DAN = 3;
n_MBON = 3;
%{
W_KD is the synaptic weight from Kenyan Cells to DANs.
The 1st row corresponds to CS+.
The 2nd row corresponds to CS-.
The first n_DAN columns correspond to attractive odors (ACV vs ETA).
The last n_DAN columns correspond to repellent o dors (OCT vs BEN).
%}
W_KD = ones(n_odor, 2*n_DAN);
%{
W_KM is the initial synaptic weight from Kenyan Cells to MBONs.
The 1st row corresponds to CS+.
The 2nd row corresponds to CS-.
%}
W_KM = ones(n_odor, n_MBON);
if is_simplified_model
    % set the weights to alpha 2 region to be 0:
    W_KD(:,[2 5]) = 0;
    W_KM(:,2) = 0;
end


fw0_wodor = eye(n_odor);
%w_pun_x = [ones(n_DAN,1), zeros(n_DAN,n_odor-1); ...
%           zeros(n_MBON,n_odor)];
w_pun_x = [[1;0;0], zeros(n_DAN,n_odor-1); ...
           zeros(n_MBON,n_odor)];
%{
W_DAN_MBON is the weight matrix between DANs and MBONs.
The rows corespond to pre-synaptic neurons;
The columns corespond to post-synaptic neurons;
%}
if is_simplified_model
    if is_keep_0_weight
        W_DAN_MBON = ... 
           [0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0;
           -1, 0,-1, 0, 0,-1;
            0, 0, 0, 0, 0, 0;
            0, 0, 1, 0, 0, 0;];
    else
        W_DAN_MBON = ... 
           [0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0;
           -1, 0,-1, 0, 0, 0;
            0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0;];
        
    end
else
    if is_keep_0_weight
        W_DAN_MBON = ... 
           [0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0;
           -1,-1,-1, 0,-1,-1;
            0, 1, 1, 0, 0, 0;
            0, 0, 1, 0, 0, 0;];
    else
        W_DAN_MBON = ... 
           [0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0;
           -1,-1,-1, 0,-1, 0;
            0, 1, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0;];
    end
end
%{
W_1 = (eye(n_DAN + n_MBON) - W_DAN_MBON')^-1;
Dx_DAN_MBON = W_1 * W_KDKM_KC';
DW_KM = ([eye(n_DAN),zeros(n_DAN, n_MBON)] * W_1 * w_pun_x)';
    
W_KDKM_KC = W_KDKM_KC + [zeros(n_odor,n_DAN),DW_KM,zeros(n_odor,n_MBON-n_DAN)];
%}
%tau_DW_KM = ones(1,n_MBON+1);
% a2 and a3 share parameters
%tau_DW_KM = ones(1,n_MBON);
%tau_DW_KM = ones(1,4);
tau_DW_KM = [1 1 1]; 
%%
mat_lu_cell = ...
    {cat(3,-inf*[W_KD(1,:),W_KM(1,:)], inf*[W_KD(1,:),W_KM(1,:)]); ...
     cat(3,-inf, 0); ...
     cat(3,-inf*w_pun_x, 0*w_pun_x); ...cat(3,-200*w_pun_x, 0*w_pun_x); ... ...
     cat(3,0*W_DAN_MBON, inf * W_DAN_MBON); ...
     cat(3,0*tau_DW_KM, inf*tau_DW_KM);%log(tau)
     cat(3,0, inf)};%log(tau)
 
for mat_i = 1:numel(mat_lu_cell)
    mat_lu_cell{mat_i}(isnan(mat_lu_cell{mat_i})) = 0 ;  
end
para_lu = parameter_vec_range(mat_lu_cell);
para_mat_cell = parameter_vec2mat(1:size(para_lu,1),mat_lu_cell);
%% set experimental conditions
[exp_para, exp_struct] = experimental_condition([]);
if size(Dx_mean,3) == 6
    % include 24hr data
    exp_para.rest_t = [3600; 3600*2-250; 3600*21-250*2];
    exp_para.session_list = exp_para.session_list([1:end,end-1,end,end-1,end]);
elseif size(Dx_mean,3) == 5
    % not include 24hr data
    exp_para.rest_t = [3600; 3600*2-250];
    exp_para.session_list = exp_para.session_list([1:end,end-1,end]);
end
[exp_para, exp_struct] = experimental_condition(exp_para);
%% pack error function
Error_weighted = @(para_vec)Error_nonlinear_activation_function_gether(...
    para_vec,mat_lu_cell,Dx_mean,Dx_SEM, ...
    @(x)Dx_steady_state_MBON_0301_2023(x,exp_struct));
Error_sq_sum = @(para_vec)nansum(reshape(Error_weighted(para_vec).^2,[],1));

is_posi_nega = para_lu==0;
Error_sq_sum_log = @(x)Error_sq_sum(para_exp_transform_simplify(x,is_posi_nega));
                
%% set initial value
if is_simplified_model
    if is_keep_0_weight
        para_0 = [-1.5;3; 24;17; 2;11;
                      3;
                      2;
                      -3;-1;-20;-20;%-0.03;-0.2;%-3;-2;-3;-1;.01;.01;
                      7;9;13; ...
                      7];%7;6;9;13;6];
    else
        para_0 = [-1.5;3; 24;17; 2;11;
                      3;
                      2;
                      -3;-2;%-0.03;-0.2;%-3;-2;-3;-1;.01;.01;
                      7;9;13; ...
                      7];%7;6;9;13;6];
    end
        %}
        %{
        para_0 = [-1.63419202817103
                2.89901356362752
                24.5888340245433
                16.3108986423069
                2.73222569298708
                12.0178779865216
                -16.8026025575586
                -4.40035246595415
                -0.0343110108689433
                -0.312319487443525
                7.43487130533623
                8.83243182555890
                12.8883569214583
                6.32003217144099]
        %}
else
    if is_keep_0_weight
        para_0 = [-1.5;-3;3; 24;18;17; 2;4;11;
                  3;
                  2;
                  -3;-2;-3;-1;-20;-20;%-0.03;-0.2;-0.2;-0.1;.01;.01;%
                  -1;-20;%0.5;.01;
                  7;10;13; ...
                  7];%7;6;9;13;6];
    else
        para_0 = [-1.5;-3;3; 24;18;17; 2;4;11;
                  3;
                  2;
                  -3;-2;-3;-1;%-0.03;-0.2;-0.2;-0.1;.01;.01;%
                  -1;%0.5;.01;
                  7;9;15; ...
                  7];%7;6;9;13;6];
    end
%}
end
for para_range_i= 3


switch para_range_i
    %{
    case 1
        options = optimset('MaxFunEvals',500000,'MaxIter',500000,...
                    'TolX',eps,'TolFun',eps, ...
                    'Display','iter','PlotFcns',@optimplotfval);

        [para_mu,fval,exitflag,output,lambda,grad,hessian] = ...
            fmincon(Error_sq_sum,...
            para_0, ...para_mu,...x,...
            [],[],[],[],para_lu(:,1),para_lu(:,2),[], ...
            options);
    %}
    case 3
        options = optimset('MaxFunEvals',500000,'MaxIter',500000,...
                    'TolX',eps,'TolFun',eps, ...
                    'Display','iter','PlotFcns',@optimplotfval);%, ...
                    %'Algorithm','active-set');

        
        
        [para_log_mu,fval,exitflag,output,lambda,grad,hessian] = ...
            fmincon(Error_sq_sum_log,...
            para_0, ...para_mu,...x,...
            [],[],[],[],-100*ones(size(para_0)),100*ones(size(para_0)),[], ...
            options);
        para_mu = para_exp_transform_simplify(para_log_mu,is_posi_nega);
        
        % Setting uplimit 100*ones(size(para_0)) to aviod the predition value to be nan.
        % Otherwise it will create 0 Square sum of error because of the nansum()

        %% confidence interval
        sample_n = 10000;
        %{
        alpha is the confidence interval level. 
        In the final manuscript, all the C.I. are at 1 sigma level. 
        %}
        %alpha = 0.05;
        alpha = normcdf(-1)*2;
        [para_log_CI, para_log_rand] = parameter_CI_sample(para_log_mu,hessian,sample_n,alpha);
        para_CI = para_exp_transform_simplify(para_log_CI,is_posi_nega);
        para_rand = para_exp_transform_simplify(para_log_rand',is_posi_nega);

        log_pdf_posterior_expx = @(x)-0.5*(Error_sq_sum_log(x)-fval);
        [fig_h,LL_mat] = plot_likelihood_ratio(para_log_mu,para_log_CI,log_pdf_posterior_expx);

    if is_use_ga
        %para_lu_ga = [para_0-5,para_0+5];
        %para_lu_ga(1:9,:) = repmat([-10 30],[9 1]);
        para_lu_ga = parameter_CI_sample(para_log_mu,hessian,sample_n,normcdf(-3));

        options = optimoptions('ga','ConstraintTolerance',1e-2,'PlotFcn', @gaplotbestf);
        [para_mu_log,fval,exitflag,output] = ...
            ga(@(x)Error_sq_sum_log(x'),...
            length(para_0), ...para_mu,...x,...
            [],[],[],[],para_lu_ga(:,1),para_lu_ga(:,2),[], ...
            options);
            para_mu = para_exp_transform_simplify(para_mu_log',is_posi_nega);
    end
    case 2
%%        %{
        para_lu_ga = [-5 5; -5 5; -5 5;...
                   0 30; 0 30; 0 30;...
                   -5 30; -5 30;-5 30; ...
                   -30 0; -30 0;
                   -3 0;-3 0;0 3;-3 0;0 3;0 3; ...
                   -3 0;-3 0;
                   5 15;5 15;5 15;5 15;];
               
        para_lu_ga = [-5 5; -5 5;...
                   0 30;  0 30;...
                   -5 30; -5 30; ...
                   -30 0; -30 0;
                   -3 0;-3 0; ...
                   5 15;5 15;5 15; ...
                   5 15;];
               
        para_lu_ga = [-50 50; -50 50;...
                   -10 30;  -10 30;...
                   -5 30; -5 30; ...
                   -30 0; -30 0;
                   -10 0; -10 0; ...
                   5 15;5 15;5 15;...
                   5 15;];
        %}
        %{
        load('/Users/flyspike/Desktop/fly model latest/figures/Dx_steady_state_nonlinear_1_05-May-2021.mat')
        para_mu_1 = para_mu;
        load('/Users/flyspike/Desktop/fly model latest/figures/Dx_steady_state_nonlinear_2_05-May-2021.mat')
        para_mu_2 = para_mu;
        para_lu_ga = (para_mu_1+para_mu_2)/2 + abs(para_mu_1-para_mu_2)*[-3 3]
        %}
        %
        %load('figures/Dx_steady_state_nonlinear_1_04-Aug-2021.mat')
        
        %para_mu_1 = para_mu;
        
        
        para_mu_1 = para_0;
        para_mu_1(abs(para_mu_1)<1e-5) = sign(para_mu_1(abs(para_mu_1)<1e-5))*0.1;
        %}
        para_lu_ga = para_mu_1*[0 1.1];
        para_lu_ga(para_mu_1<0,:) = para_lu_ga(para_mu_1<0,[2 1])
        %para_lu_ga(abs(para_lu_ga(:,2))<1e-5,2) = 0.1
        %{
        para_lu_ga = ...
        [-3.23447180675298,-0.527917852221858;-8.79783764168201,11.2308475569653;-19.6083117803271,19.2515346274710;
         26.0091627470069,28.2985339578084;-0.733742936769058,33.8015833045811;13.6769350008379,21.1746053701930;
         2.03350248762839,3.74070514954688;-2.11249950763431,20.7054140899136;-6.91426397724747,27.1103470917760;
         -40,0;-10,-4;
         -0.0465439869466340,-0.0285884182747124;-0.616866210373795,0;0,0.145192259034053;
         -0.503928203035706,-0.0150512835459244;0,1.53124996592269;0,0.5;
         -0.944807249260022,0;-0.5,0;
         6.92035377043877,7.89774898969157;4.76156679359778,12.07558287064444;8.48289079832680,9.49433641116539;7.68930381435289,15;5.25226074168680,7.15911496059413]
        %}
        options = optimoptions('ga','ConstraintTolerance',1e-2,'PlotFcn', @gaplotbestf);
        [para_mu,fval,exitflag,output] = ...
    ga(Error_sq_sum,...
    length(para_0), ...para_mu,...x,...
    [],[],[],[],para_lu_ga(:,1),para_lu_ga(:,2),[], ...
    options);
    para_mu = para_mu';

    %{
    para_mu = [-1.56921560560083
2.71334534078927
24.3871857946204
16.3681971141920
2.60075346034134
11.8322174361735
-15.6061016380302
-4.10054185609696
-0.0306448850580688
-0.308715060643284
7.41312524293108
8.84246646297281
12.9653089403717
6.20970844077093]
    %}
end

    %%
    fig_h = [0 0];
    [error_weighted,Dx_DAN_MBON] = Error_weighted(para_mu);

    for fig_i = 1:2
        fig_h(fig_i) = figure;
        axis_h = zeros(3,2);
        for axis_i = 1:numel(axis_h)
            axis_h(axis_i)=subplot(size(axis_h,2),size(axis_h,1),axis_i); hold on;
        end
        axis_h=axis_h';

        %para_mu_curr = para_mu;
        %if fig_i == 2
        %    para_mu_curr(1:3,1) = para_mu(7:9,1);
        %end
        %Dx_fun = @(para)Dx_steady_state_MBON_ReLu_detail(para,exp_struct);
        plot_Dx_data_new(axis_h,Dx_mean(:,:,:,fig_i),Dx_SEM(:,:,:,fig_i),Dx_DAN_MBON(:,:,:,fig_i));
    end
Error_sq_sum(para_mu)
%ans = 160.3898
%%
%{
para_mat_cell = [parameter_vec2mat(para_mu,mat_lu_cell);
                 parameter_vec2mat(para_CI(:,1),mat_lu_cell);
                 parameter_vec2mat(para_CI(:,2),mat_lu_cell)];
%}
para_mat_cell = parameter_vec2mat(para_mu,mat_lu_cell);

fig_name = fullfile('figures',['Dx_steady_state_nonlinear_',num2str(para_range_i),...
    '_',date,'_',num2str(3-is_simplified_model),'modules']);
save(fig_name,'para_mu','para_CI','mat_lu_cell','para_rand')
%save(fig_name,'para_mu','mat_lu_cell')
saveas(fig_h(1),[fig_name,'_1'])
saveas(fig_h(2),[fig_name,'_2'])
print(fig_h(1),[fig_name,'_1'],'-dtiff','-r300')
print(fig_h(2),[fig_name,'_2'],'-dtiff','-r300')
%{
xlswrite(fullfile('figures','para_list.xlsx'), ...
         para_xls_cell(para_mat_cell), ...
         ['Dx_SS_NL_',num2str(para_range_i),'_',date])
%}
if is_simplified_model
    xls_start = 'C3';
else
    xls_start = 'F3';
end
%%
xlswrite(fullfile('figures','para_compare.xlsx'), ...
         para_mu_CI_xls(para_mu,para_CI,mat_lu_cell), ...
         ['Dx_SS_NL_',num2str(para_range_i),'_',date], ...
         xls_start)

end
%%
% data_point_n = sum(reshape(~isnan(Dx_mean),[],1));
% fun_nonnan = @(x)reshape(x(~isnan(x)),[],1);
% modelfun = @(b,x)fun_nonnan(Error_weighted(para_mu.*abs(b)));
%      mdl = fitnlm(1:data_point_n,zeros(1,data_point_n), ...
%                   modelfun,ones(size(para_mu)));
% beta = table2array(mdl.Coefficients(:,1));
end

%%
%{
rand_n = 1000;
rand_Error = zeros(1,rand_n);
rand_para_mu = zeros(size(para_mu,1),rand_n);

for rand_i = 1:rand_n
    para_mu_curr = para_mu.*(1+randn(size(para_mu))); 
    rand_Error(rand_i) = Error_sq_sum(para_mu_curr);
    rand_para_mu(:,rand_i) = para_mu_curr;
end

figure,
hist(rand_Error)
%}