function para_mat_cell = reset_parameter_by_valance(para_mu,mat_lu_cell,odor_valence)
para_mat_cell = parameter_vec2mat(para_mu,mat_lu_cell);
%{
W_KDKM_KC = [para_mat_cell{1}(1:6);para_mat_cell{1}(1:6)];%
W_DAN_MBON = para_mat_cell{4};

n_DAN = 3;
n_MBON = size(W_DAN_MBON,1) - n_DAN;
W_1 = (eye(n_DAN + n_MBON) - W_DAN_MBON')^-1;
Firing_Temp = W_1 * W_KDKM_KC';
%}
KC_adapt = 0.05;
Firing_Temp = zeros(6,1);
Firing_Temp(1:3,:) = odor_valence;
Firing_Temp(4,:) = 31.6;
Firing_Temp(5,:) = 5.8;
Firing_Temp(6,:) = 17.8;

W_DAN_MBON = para_mat_cell{4};
W_KDKM_KC = (eye(size(W_DAN_MBON))-W_DAN_MBON')*Firing_Temp*(2/(1+exp(-KC_adapt*5)));
para_mat_cell{1} = W_KDKM_KC';

end