clc,clear
addpath(fullfile('..','..','data_and_parameters'))
addpath(fullfile('..','model_functions'))
addpath(fullfile('..','output_functions'))
addpath('customcolormap')

file_parameter_list = ...
    {'Dx_steady_state_nonlinear_3_27-Mar-2023_2modules.mat'  
     'Dx_steady_state_nonlinear_3_27-Mar-2023_3modules.mat'};

figure_type_list = {'ACV_EtA','OCT_BEN'};
for figure_type_i= 1:length(figure_type_list)
for file_i = 1:2
    fun_plotPanel_5j(fullfile('..','..','data_and_parameters',file_parameter_list{file_i}), ...
                     fullfile('figures',['figure_5j_',figure_type_list{figure_type_i},file_parameter_list{file_i}(29:end-4)]), ...
                     figure_type_list{figure_type_i});
end
end