function [error_weighted,Dx_DAN_MBON] = Error_nonlinear_activation_function_gether(para_vec,mat_lu_cell,Dx_mean,Dx_SEM,Dx_steady_state_nonlinear)
para_mat_cell_all = parameter_vec2mat(para_vec,mat_lu_cell);
%% ACV vs ETA
para_mat_cell = para_mat_cell_all;
para_mat_cell{1} = para_mat_cell{1}(:,1:6);

Dx_DAN_MBON_ACV_ETA = Dx_steady_state_nonlinear(para_mat_cell);
%% OCT vs BEN
para_mat_cell = para_mat_cell_all;
para_mat_cell{1} = para_mat_cell{1}(:,[7:9 4:6]);

Dx_DAN_MBON_OCT_BEN = Dx_steady_state_nonlinear(para_mat_cell);
%%
Dx_DAN_MBON = cat(4,Dx_DAN_MBON_ACV_ETA, Dx_DAN_MBON_OCT_BEN);

error_weighted = (Dx_DAN_MBON-Dx_mean)./Dx_SEM;
%ES = nansum(reshape(error_weighted.^2,[],1));
end