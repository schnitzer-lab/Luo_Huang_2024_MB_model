%clc,
clear
addpath(fulfile('..','model_functions'))
addpath(fulfile('..','output_functions'))
%%
% if is_simplified_model == true : 2-module model
% if is_simplified_model == false: 3-module model
for is_simplified_model = [true false]
    % some of the synaptic weights is close to 0. 
    % We can try to set them to be 0 by let is_keep_0_weight = false;
    is_keep_0_weight = true;%false;
    % The parameter is_use_ga set whether we use the genetic algorithm to
    % confirm the parameter is not a local minimum.
    is_use_ga = false;% true
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
    w_pun_x = [[1;0;0], zeros(n_DAN,n_odor-1); ...
               zeros(n_MBON,n_odor)];
    %{
    W_DAN_MBON is the weight matrix between DANs and MBONs.
    The rows corespond to pre-synaptic neurons;
    The columns corespond to post-synaptic neurons;
    %}
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
    if is_simplified_model
        W_DAN_MBON([2 5],:) = 0;
        W_DAN_MBON(:,[2 5]) = 0;
    end
    % tau_DM_KM are the Timeconstants:
    % τMBON-STM (KC to MBON-γ1pedc>α/β),
    % τMBON-STM (KC to MBON-α'2α2 and KC to MBON-α3)(t ≤ 3 hrs)
    % τMBON-LTM (KC to MBON-α'2α2 and KC to MBON-α3)(t > 3 hrs)
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
                          -3;-1;-20;-20;
                          7;9;13; ...
                          7];
        else
            para_0 = [-1.5;3; 24;17; 2;11; ...
                          3;
                          2;
                          -3;-2;
                          7;9;13; ...
                          7];
        end
            
    else
        if is_keep_0_weight
            para_0 = [-1.5;-3;3; 24;18;17; 2;4;11;
                      3;
                      2;
                      -3;-2;-3;-1;-20;-20;
                      -1;-20;
                      7;10;13; ...
                      7];
        else
            para_0 = [-1.5;-3;3; 24;18;17; 2;4;11; ...
                      3;
                      2;
                      -3;-2;-3;-1;
                      -1;
                      7;9;15; ...
                      7];
        end
    end
    %%
    options = optimset('MaxFunEvals',500000,'MaxIter',500000,...
                'TolX',eps,'TolFun',eps, ...
                'Display','iter','PlotFcns',@optimplotfval);

    [para_log_mu,fval,exitflag,output,lambda,grad,hessian] = ...
        fmincon(Error_sq_sum_log,...
        para_0, ...
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
    %% Use genetic algorithm to scan larger parameter space.
    % This step avoids the parameters fall to local minimum.
    if is_use_ga
        para_lu_ga = parameter_CI_sample(para_log_mu,hessian,sample_n,normcdf(-6)*2);

        options = optimoptions('ga','ConstraintTolerance',1e-2,'PlotFcn', @gaplotbestf);
        [para_mu_log,fval,exitflag,output] = ...
            ga(@(x)Error_sq_sum_log(x'),...
            length(para_0), ...para_mu,...x,...
            [],[],[],[],para_lu_ga(:,1),para_lu_ga(:,2),[], ...
            options);
            para_mu = para_exp_transform_simplify(para_mu_log',is_posi_nega);
    end
    

    %% plot data
    fig_h = [0 0];
    [error_weighted,Dx_DAN_MBON] = Error_weighted(para_mu);

    for fig_i = 1:2
        fig_h(fig_i) = figure;
        axis_h = zeros(3,2);
        for axis_i = 1:numel(axis_h)
            axis_h(axis_i)=subplot(size(axis_h,2),size(axis_h,1),axis_i); hold on;
        end
        axis_h=axis_h';
        plot_Dx_data_new(axis_h,Dx_mean(:,:,:,fig_i),Dx_SEM(:,:,:,fig_i),Dx_DAN_MBON(:,:,:,fig_i));
    end
    if is_simplified_model
        delete(axis_h(:,2))
    end

    Error_sq_sum(para_mu)
    %% save figures
    para_mat_cell = parameter_vec2mat(para_mu,mat_lu_cell);
    
    fig_name = fullfile('figures',['Dx_steady_state_nonlinear_3',...
        '_',date,'_',num2str(3-is_simplified_model),'modules']);
    save(fig_name,'para_mu','para_CI','mat_lu_cell','para_rand')

    saveas(fig_h(1),[fig_name,'_1'])
    saveas(fig_h(2),[fig_name,'_2'])
    print(fig_h(1),[fig_name,'_1'],'-dtiff','-r300')
    print(fig_h(2),[fig_name,'_2'],'-dtiff','-r300')
    %% write parameters into xls file. 
    if is_simplified_model
        xls_start = 'F3';
    else
        xls_start = 'C3';
    end
    xlswrite(fullfile('figures','para_compare.xlsx'), ...
             para_mu_CI_xls(para_mu,para_CI,mat_lu_cell), ...
             ['Dx_SS_NL_3','_',date], ...
             xls_start)
end


