function [Dx_DAN_MBON,Dx_0] = solve_DAN_MBON(Dx_KC, x_punish,W_KM, w_punish, W_mat, fun_act, Ab_DAN_MBON)
%{
Dx_KC: KC input, n_odor*1 column vector
x_punish: punishment input, single value
W_KM: weight from KCs to MBON, n_odor*n_MBON matrix
w_punish: weight of punishment n_DAN*1 column vector
W_mat = [0 0; W_MD W_MM]: weight matrix between MBON and DAN 
        (n_DAN+n_MBON)*(n_DAN+n_MBON) matrix
fun_act: activation function
Ab_DAN_MBON: The amplitue and bias of DAN and MBON. (n_DAN+n_MBON)*2 matrix
%}

B_DAN_MBON = Ab_DAN_MBON(:,1).*fun_act(Ab_DAN_MBON(:,2));
WDx_KC_pun = W_KM'*Dx_KC + w_punish*x_punish;
WDx_KC_pun_b = WDx_KC_pun + Ab_DAN_MBON(:,2);
fun_Dx = @(Dx)Ab_DAN_MBON(:,1).*fun_act(WDx_KC_pun_b + W_mat'*Dx) - B_DAN_MBON - Dx;

Dx_0 = (eye(size(W_mat))-W_mat)^-1 * WDx_KC_pun;
for re_i = 1:10
Dx_0 = fun_Dx(Dx_0) + Dx_0;
end
Dx_DAN_MBON  = Dx_0;
%Dx_DAN_MBON  = fsolve(fun_Dx,Dx_0);
end

