clc,clear all
addpath(fullfile('..','..','data_and_parameters'))
addpath(fullfile('..','model_functions'))
addpath(fullfile('..','output_functions'))
%Dx_steady_state_MBON_revision = @(x,y)Dx_steady_state_MBON_0301_2023(x,y);
%%
%{ 
This script plot figure 5c (3-module model) and S10j (2-module model).

% if is_simplified_model == true : 2-module model
% if is_simplified_model == false: 3-module model
%}
for is_simplified_model = [true false]
    %% load original data
    data_ori_filename = fullfile('..','..','data_and_parameters','Imaging_24hr_data');
    [Dx_mean, Dx_SEM] = load_original_data(data_ori_filename, is_simplified_model);
    %{
    Dimension 1: PPL1-v1pedc, PPL1-a'2a2, PPL1-a3, MBON-v1pedc, MBON-a2sc, MBON-a3
    Dimension 2: CS+, CS-
    Dimension 3: Pre-traing, after 3x training, after 6x training, 1hr memory, 3hr, 24hr
    Dimension 4: ACV vs ETA, OCT vs BEN
    %}
    
    %%
    % create the default experimental condition list
    [exp_para, exp_struct] = experimental_condition([]);
    exp_para.session(3).n = 3; % training cycle
    [exp_para, exp_struct] = experimental_condition(exp_para);
    
    exp_para.rest_t = [3600-250-300; 3600*2-250; 3600*21-250];
    exp_para.session_list = exp_para.session_list([1:end,end-1,end,end-1,end]);
    [exp_para, exp_struct] = experimental_condition(exp_para);
    exp_struct(4).t_length = 300;
    exp_struct(16).t_length = 300;
    exp_struct(20).t_length = 300;
    exp_struct(32).t_length = 300;

%% load parameters
    if is_simplified_model
        load('Dx_steady_state_nonlinear_3_27-Mar-2023_2modules') %%Extended Data Fig10j
    else
        load('Dx_steady_state_nonlinear_3_27-Mar-2023_3modules') %% Fig5c
    end

    Dx_DAN_MBON = zeros([size(Dx_mean), size(para_rand,2)]);
    %{
    Dimension 1: PPL1-v1pedc, PPL1-a'2a2, PPL1-a3, MBON-v1pedc, MBON-a2sc, MBON-a3
    Dimension 2: CS+, CS-
    Dimension 3: Pre-traing, after 3x training, after 6x training, 1hr memory, 3hr, 24hr
    Dimension 4: ACV vs ETA, OCT vs BEN
    Dimention 5: parameter sample i
    %}
    for para_i = 1:size(para_rand,2)
        
        para_mat_cell = parameter_vec2mat(para_rand(:,para_i),mat_lu_cell);
        para_mat_cell_ACV_ETA = para_mat_cell;
        para_mat_cell_ACV_ETA{1} = para_mat_cell{1}(1,1:6);
    
        para_mat_cell_OCT_BEN = para_mat_cell;
        para_mat_cell_OCT_BEN{1} = para_mat_cell{1}(1,[7:9 4:6]);
        
        Dx_DAN_MBON(:,:,:,:,para_i) = cat(4, ...
            Dx_steady_state_MBON_0301_2023(para_mat_cell_ACV_ETA,exp_struct), ...
            Dx_steady_state_MBON_0301_2023(para_mat_cell_OCT_BEN,exp_struct));
        %{
        [~, Dx_DAN_MBON(:,:,:,:,para_i)] = ...
            Error_nonlinear_activation_function_gether(...
            para_rand(:,para_i),mat_lu_cell,Dx_mean,Dx_SEM, ...
            @(x)Dx_steady_state_MBON_0301_2023(x,exp_struct));
        %}
    end
    Dx_DAN_MBON_median_eb = prctile(Dx_DAN_MBON,normcdf([0 -1 1])*100,5);

%%
    [fig_h, error_bar_h, fit_h] = plot_Dx_data(Dx_mean,Dx_SEM,Dx_DAN_MBON_median_eb);
    
end
