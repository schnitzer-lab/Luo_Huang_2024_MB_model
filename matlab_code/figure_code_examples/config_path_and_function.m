addpath(fullfile('..','..','data_and_parameters'))
addpath(fullfile('..','model_functions'))
addpath(fullfile('..','output_functions'))
Dx_steady_state_MBON_revision = @(x,y)Dx_steady_state_MBON_0301_2023(x,y);