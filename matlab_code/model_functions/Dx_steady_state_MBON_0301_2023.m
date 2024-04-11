function [Dx_DAN_MBON,data_all] = Dx_steady_state_MBON_0301_2023(para_mat_cell,exp_struct)
%{
This function simulate the learning process by recurrent model.

para_mat_cell: the parameters in cell format.
exp_struct: the struct of training protocol in experiments.
%}

%% parameters from experiments and literature
w_pun = [27.85;0;11.38;0;0;0];
% The firing rate of PPL1-gamma1 neuron increased by 27.85 +/- 3.28 Hz during shock, 
% while PPL1-a3 increased by 11.38  +/- 2.53 Hz.

B_MBON = [35.2; 9.0; 11.2];% new baseline value % from file OdorData(revision).xlsx
max_MBON = B_MBON +[36.46; 8.9; 19.96];% ~ mean + std

KC_adapt = 0.05;% the adaptation 
%% input parameters

fw0 = para_mat_cell{2};%*eye(2);%
% During the imaging, there was no punishment input. The punishment only
% affect the training. So we can wrap the fwDt_w_pun_odor.
% we also need to notice here used the linear approximation on DANs.

% fwDt_w_pun = para_mat_cell{3};%%%%%%
fwDt = para_mat_cell{3}(1);%%%%%%
fwDt_w_pun = fwDt * w_pun;

W_KDKM_KC = [para_mat_cell{1};para_mat_cell{1}];%

W_DAN_MBON = para_mat_cell{4};
%tau_DW_KM = exp(para_mat_cell{5});
%tau_odor_adaptation = exp(para_mat_cell{6});
tau_DW_KM = para_mat_cell{5};
tau_odor_adaptation = para_mat_cell{6};

n_odor = size(W_KDKM_KC,1);
n_DAN = 3;
n_MBON = size(W_DAN_MBON,1) - n_DAN;
%{
training_n = training_para.training_n;%3;%
training_interval = training_para.training_interval;%135;%
t_between_tests = (5+120)*2+(30+training_interval)*2*training_n;%
t_rest = training_para.t_rest;%3600;%
%}
% training cycle number

%%

Ab_DAN_MBON = [ones(size(W_DAN_MBON,1),1), ...
               [0; 0; 0;B_MBON]];
           
%fun_act = @(x)[x(1:n_DAN,:);max(x(n_DAN+1:end,:),0)];
% the activation function of DANs  is x, 
% the activation function of MBONs is piece-wise linear, 
fun_act = @(x)[x(1:n_DAN,:);min(max(x(n_DAN+1:end,:),0),max_MBON)];
%}
%%
W_KDKM_KC_0 = W_KDKM_KC;
DW_KM = 0;
w_odor_start = [1;1];
imaging_ind = 1;
Dx_DAN_MBON = zeros(size(W_DAN_MBON,1),...
                    length(exp_struct(1).odor),...
                    sum([exp_struct.imaging_i]==1));
if nargout == 2
    data_all = struct( ...
        'w_odor_mean',zeros([size(w_odor_start), length(exp_struct)]), ...
        'w_odor_end' ,zeros([size(w_odor_start), length(exp_struct)]), ...
        'W_KDKM_KC'  ,zeros([size(   W_KDKM_KC), length(exp_struct)]), ...
        'Dx_KC_mean' ,zeros([size(w_odor_start), length(exp_struct)]), ...
        'Dx_DAN_MBON',zeros([size(Dx_DAN_MBON,1),1,length(exp_struct)]));
end


% a2 and a3 share parameters
% tau_DW_KM_curr = tau_DW_KM([1 2 2]);
last_training_i = find(strcmp({exp_struct.name},'training'), 1, 'last' );
for exp_i = 1:length(exp_struct)
    
    w_odor_end = w_odor_start.*exp(-KC_adapt*exp_struct(exp_i).t_length*exp_struct(exp_i).odor);

    w_odor_end = [1;1] - ...
        ([1;1]-w_odor_end)*exp(-exp_struct(exp_i).t_length/tau_odor_adaptation);
    w_odor_mean = (w_odor_start+w_odor_end)/2;%[1;1];%
    w_odor_start = w_odor_end;
    
    Dx_KC_mean = w_odor_mean.*exp_struct(exp_i).odor;
    Dx_DAN_MBON_curr = solve_DAN_MBON( ...
        Dx_KC_mean, exp_struct(exp_i).punishment,...diag(Dx_KC_mean)
        W_KDKM_KC, w_pun, ...
        W_DAN_MBON, fun_act, Ab_DAN_MBON);

    %{
    We considered 3 adjacent 30-s training bouts as the standard training condition (Fig. 3d). 
    The A_AH(0) of imaging bout and extinction bouts are normalized to them.
    If the time length of other bouts are not 30 s, the A_AH(0) should be
    changed propotionally. 
    %}
    fw0_curr = fw0*exp_struct(exp_i).t_length/(3*30);
    fwDt_w_pun_curr = fwDt_w_pun(1:n_DAN,1)*exp_struct(exp_i).t_length/(3*30); %fwDt_w_pun(1:n_DAN,:)
    
    DW_KM = DW_KM + ((fw0_curr*(W_KDKM_KC(:,1:n_DAN)' * Dx_KC_mean + ...diag(Dx_KC_mean)
           W_DAN_MBON(n_DAN+1:end,1:n_DAN)'*Dx_DAN_MBON_curr(n_DAN+1:end,:)) ...
         + fwDt_w_pun_curr * exp_struct(exp_i).punishment)*Dx_KC_mean')';%diag(Dx_KC_mean)
    
    %{%exp law
    %{
    t_shift_tau = 3600*3;
    
    if exp_struct(exp_i).t_length <= t_shift_tau
        % consider long-term consolidation
        % a2 and a3 share parameters
        tau_DW_KM_curr = tau_DW_KM([1 2 2]);
        DW_KM = DW_KM.*exp(-exp_struct(exp_i).t_length./tau_DW_KM_curr);
    else
        %tau_DW_KM_curr = tau_DW_KM([1 2 2]);
        %DW_KM = DW_KM.*exp(-t_shift_tau./tau_DW_KM_curr);
        tau_DW_KM_curr = tau_DW_KM([1 3 3]);
        %DW_KM = DW_KM.*exp(-(exp_struct(exp_i).t_length-t_shift_tau)./tau_DW_KM_curr);
        DW_KM = DW_KM.*exp(-exp_struct(exp_i).t_length./tau_DW_KM_curr);
    end
    %}
    t_shift_tau = 3600*3;
    t_cum = sum([exp_struct(last_training_i+1:exp_i).t_length]);
    if t_cum <= t_shift_tau
        % consider long-term consolidation
        % a2 and a3 share parameters
        tau_DW_KM_curr = tau_DW_KM([1 2 2]);
        DW_KM = DW_KM.*exp(-exp_struct(exp_i).t_length./tau_DW_KM_curr);
    elseif t_shift_tau > t_cum - exp_struct(exp_i).t_length
        t_shift_rest = t_shift_tau - (t_cum - exp_struct(exp_i).t_length);
        tau_DW_KM_curr = tau_DW_KM([1 2 2]);
        DW_KM = DW_KM.*exp(-t_shift_rest./tau_DW_KM_curr);
        tau_DW_KM_curr = tau_DW_KM([1 3 3]);
        DW_KM = DW_KM.*exp(-(exp_struct(exp_i).t_length-t_shift_rest)./tau_DW_KM_curr);
        %DW_KM = DW_KM.*exp(-exp_struct(exp_i).t_length./tau_DW_KM_curr);
    else
        tau_DW_KM_curr = tau_DW_KM([1 3 3]);
        DW_KM = DW_KM.*exp(-exp_struct(exp_i).t_length./tau_DW_KM_curr);
    end


    %}
    %{
    % power law     
	DW_KM(:,1) = DW_KM(:,1)*(1+exp_struct(exp_i).t_length/tau_DW_KM(1))^(-tau_DW_KM(3));
    DW_KM(:,[2 3]) = DW_KM(:,[2 3])*(1+exp_struct(exp_i).t_length/tau_DW_KM(2))^(-tau_DW_KM(4));
    %}
    W_KDKM_KC = W_KDKM_KC_0 + [zeros(n_odor,n_DAN),DW_KM,zeros(n_odor,n_MBON-n_DAN)];
    % ReLu on Weights from KC to MBON
    % W_KDKM_KC(:,4:end) = max(W_KDKM_KC(:,4:end),0);
    % softplus on Weights from KC to MBON
    % W_KDKM_KC(:,4:end) = log(1+exp(W_KDKM_KC(:,4:end)));
    % consider the parallel inhibition effects and DPM neuron, the overall weight
    % from KCs to MBONs can be negative.
    
    if exp_struct(exp_i).imaging_i > 0
        Dx_DAN_MBON(:,exp_struct(exp_i).imaging_i,imaging_ind) = Dx_DAN_MBON_curr;
        imaging_ind = imaging_ind + (exp_struct(exp_i).imaging_i==2);
    end
    %%
    if nargout == 2
        data_all.w_odor_mean(:,:,exp_i) = w_odor_mean;
        data_all.w_odor_end (:,:,exp_i) = w_odor_end;
        data_all.W_KDKM_KC  (:,:,exp_i) = W_KDKM_KC;
        data_all.Dx_KC_mean (:,:,exp_i) = Dx_KC_mean;
        data_all.Dx_DAN_MBON(:,:,exp_i) = Dx_DAN_MBON_curr;
    end
    
end

end


function [Dx_DAN_MBON,Dx_0] = solve_DAN_MBON(Dx_KC, x_punish,W_KM, w_punish, W_mat, fun_act, Ab_DAN_MBON)
%{
This function solve the nonlinear function of MBON

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