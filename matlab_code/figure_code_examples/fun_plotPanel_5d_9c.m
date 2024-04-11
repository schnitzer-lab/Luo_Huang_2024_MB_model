function fun_plotPanel_5d_9c(file_parameter,figure_name,figure_type)


%load('/Users/flyspike/Desktop/Fly leanring paper/fly model with CI/version2 - a2 a3 equal tau/figures/Dx_steady_state_nonlinear_3_18-Aug-2021')
%load('figures/Dx_steady_state_nonlinear_3_18-Aug-2021')
load(file_parameter)
addpath('customcolormap\')

%-----loading model fitted parameters
%para_mat_cell = parameter_vec2mat(x_lu,mat_lu_cell);


%{
fw0 = para_mat_cell{2};%*eye(2);%
fwDt = para_mat_cell{3}(1);%%%%%%

%w_pun = [26.81;0;10.00;0;0;0];
w_pun = [28.52;0;12.18;0;0;0];
% The firing rate of PPL1-gamma1 neuron increased by 28.52 +/- 3.28 Hz during shock, 
% while PPL1-a3 increased by 12.18  +/- 2.57 Hz.

para_mat_cell{3}(:,1) = fwDt * w_pun;

W_KDKM_KC = [para_mat_cell{1}(1:6);para_mat_cell{1}(1:6)];%
W_DAN_MBON = para_mat_cell{4};
tau_DW_KM = exp(para_mat_cell{5});
tau_odor_adaptation = exp(para_mat_cell{6});

n_odor = size(W_KDKM_KC,1);
n_DAN = 3;
n_MBON = size(W_DAN_MBON,1) - n_DAN;
W_1 = (eye(n_DAN + n_MBON) - W_DAN_MBON')^-1;
Firing_Temp = W_1 * W_KDKM_KC';

%% ----------ACV-ETA
KC_adapt = 0.05;
Firing_Temp(4,:) = 31.6*(2/(1+exp(-KC_adapt*5)));
Firing_Temp(5,:) = 5.8*(2/(1+exp(-KC_adapt*5)));
Firing_Temp(6,:) = 17.8*(2/(1+exp(-KC_adapt*5)));

%% -----Modified connection matrix
% W_DAN_MBON(4,1:6) = 0; % disconnect MBON-v1 feedbacks
% W_DAN_MBON(5,1:6) = 0; % disconnect MBON-a2sc feedbacks
% W_DAN_MBON(6,1:6) = 0; % disconnect MBON-a3 feedbacks

W_1 = (eye(n_DAN + n_MBON) - W_DAN_MBON')^-1;
W_KDKM_KC = inv(W_1)*Firing_Temp;
W_KDKM_KC = W_KDKM_KC';

%-----Training paradigm parameters
shock_strength = 1;
fwDt_w_pun = fwDt * w_pun * shock_strength;
%}
% create the defualt experimental condition listexp_para
[exp_para, exp_struct] = experimental_condition([]);

exp_para.session_list = ...
    {'imaging';
    'training';
    'rest';
    'imaging'};

switch figure_type
    case '5d'
        exp_para.session(3).t_ISI = 135; % interval
        exp_para.session(3).n = 1; % training cycle % Figure 5d
        odor_valence_list = -4:1:4;
        rest_t_list = (300:300:3600)-135;
    case 'S9c'
        exp_para.session(3).t_ISI = 135; % interval
        exp_para.session(3).n = 6; % training cycle % S-figure 8c
        odor_valence_list = -8:1:8;
        rest_t_list = (300:300:86400)-135;
    case '2train360'
        exp_para.session(3).t_ISI = 360; % interval
        exp_para.session(3).n = 2; 
        odor_valence_list = -8:1:8;
        rest_t_list = (300:300:86400)-135;
    case '2train3600'
        exp_para.session(3).t_ISI = 3600; % interval
        exp_para.session(3).n = 2; 
        odor_valence_list = -8:1:8;
        rest_t_list = (300:300:86400)-135;
    otherwise
        error('figure_type should be 5d, S9c, 2train360 or 2train3600.')
end
[exp_para, exp_struct] = experimental_condition(exp_para);

%%
Dx_DAN_MBON_Matrix_CSp = cell(6,1);
Dx_DAN_MBON_Matrix_CSm = cell(6,1);
Dx_DAN_MBON_Matrix_delta = cell(6,1);
MemoryTrace_CSp = cell(6,1);
MemoryTrace_CSm = cell(6,1);
MemoryTrace_delta = cell(6,1);

indx = 1;
%W_KDKM_KC_0 = W_KDKM_KC;

for odor_valence = odor_valence_list%-4:1:4 % Figure 5d

% for odor_valence = -8:1:8 % S-figure 8c
%{
    NoValence = W_1*W_KDKM_KC_0';
    NoValence(1:3,:) = 1 * odor_valence*(2/(1+exp(-KC_adapt*5)));
    W_KDKM_KC_0 = (inv(W_1)*NoValence)';
%}
    para_mat_cell = reset_parameter_by_valance(para_mu,mat_lu_cell,odor_valence);
    indxx = 1;
    
    for rest_t = rest_t_list%(300:300:3600)-135 % Figure 5d
%     for rest_t = (300:300:86400)-135 % S-figure 8c

        %     for rest_t = 10250
        
        exp_para.rest_t = rest_t;
        
        
        [exp_para, exp_struct] = experimental_condition(exp_para);
        exp_struct(4).t_length = 300;
        %{

        B_MBON = [35.2; 9.0; 11.2];% new baseline value
        
        Ab_DAN_MBON = [ones(size(W_DAN_MBON,1),1), ...
            [0; 0; 0;B_MBON]];
        
        max_MBON = B_MBON +[36.46; 8.9; 19.96];% ~ mean + std
        
        
        fun_act = @(x)[x(1:n_DAN,:);min(max(x(n_DAN+1:end,:),0),max_MBON)];
        
                
        W_KDKM_KC = W_KDKM_KC_0;
        DW_KM = 0;
        w_odor_start = [1;1];
        imaging_ind = 1;
        Dx_DAN_MBON = zeros(size(W_DAN_MBON,1),...
            length(exp_struct(1).odor),...
            sum([exp_struct.imaging_i]==1));
        
        data_all = struct( ...
            'w_odor_mean',zeros([size(w_odor_start), length(exp_struct)]), ...
            'w_odor_end' ,zeros([size(w_odor_start), length(exp_struct)]), ...
            'W_KDKM_KC'  ,zeros([size(W_KDKM_KC), length(exp_struct)]), ...
            'Dx_KC_mean' ,zeros([size(w_odor_start), length(exp_struct)]), ...
            'Dx_DAN_MBON',zeros([size(Dx_DAN_MBON,1),1,length(exp_struct)]));
        
        tau_DW_KM_curr = tau_DW_KM([1 2 2]);
        for exp_i = 1:length(exp_struct)
            KC_adapt = 0.05;
            w_odor_end = w_odor_start.*exp(-KC_adapt*exp_struct(exp_i).t_length*exp_struct(exp_i).odor);
            w_odor_end = [1;1] - ...
                ([1;1]-w_odor_end)*exp(-1*exp_struct(exp_i).t_length/tau_odor_adaptation);
            w_odor_mean = (w_odor_start+w_odor_end)/2;%[1;1];%
            w_odor_start = w_odor_end;
            
            Dx_KC_mean = w_odor_mean.*exp_struct(exp_i).odor;
            Dx_DAN_MBON_curr = solve_DAN_MBON( ...
                Dx_KC_mean, exp_struct(exp_i).punishment,...diag(Dx_KC_mean)
                W_KDKM_KC, w_pun, ...
                W_DAN_MBON, fun_act, Ab_DAN_MBON);
            
            fw0_curr = fw0*exp_struct(exp_i).t_length/(3*30);
            fwDt_w_pun_curr = fwDt_w_pun(1:n_DAN,1)*exp_struct(exp_i).t_length/(3*30); %fwDt_w_pun(1:n_DAN,:)
            DW_KM = DW_KM + ((fw0_curr*(W_KDKM_KC(:,1:n_DAN)' * Dx_KC_mean + ...diag(Dx_KC_mean)
                W_DAN_MBON(n_DAN+1:end,1:n_DAN)'*Dx_DAN_MBON_curr(n_DAN+1:end,:)) ...
                + fwDt_w_pun_curr * exp_struct(exp_i).punishment)*Dx_KC_mean')';%diag(Dx_KC_mean)
            
            if exp_struct(exp_i).t_length > 3600*3
                tau_DW_KM_curr1 = tau_DW_KM([1 3 3]);
                % consider long-term consolidation
                DW_KM1 = DW_KM.*exp(-3600*3./tau_DW_KM_curr);
                DW_KM = DW_KM1.*exp(-(exp_struct(exp_i).t_length-3600*3)./tau_DW_KM_curr1);
            else
                DW_KM = DW_KM.*exp(-exp_struct(exp_i).t_length./tau_DW_KM_curr);
            end
            
            W_KDKM_KC = W_KDKM_KC_0 + [zeros(n_odor,n_DAN),DW_KM,zeros(n_odor,n_MBON-n_DAN)];
            
            if exp_struct(exp_i).imaging_i > 0
                Dx_DAN_MBON(:,exp_struct(exp_i).imaging_i,imaging_ind) = Dx_DAN_MBON_curr;
                imaging_ind = imaging_ind + (exp_struct(exp_i).imaging_i==2);
            end
            %%
            data_all.w_odor_mean(:,:,exp_i) = w_odor_mean;
            data_all.w_odor_end (:,:,exp_i) = w_odor_end;
            data_all.W_KDKM_KC  (:,:,exp_i) = W_KDKM_KC;
            data_all.Dx_KC_mean (:,:,exp_i) = Dx_KC_mean;
            data_all.Dx_DAN_MBON(:,:,exp_i) = Dx_DAN_MBON_curr;
            
        end
        %}
        [Dx_DAN_MBON,data_all] = Dx_steady_state_MBON_0301_2023(para_mat_cell,exp_struct);
        for neuron_inx = 1:6
            for ind = 1:size(Dx_DAN_MBON,3)
                Dx_DAN_MBON_Matrix_CSp{neuron_inx,1}(indxx,ind) = Dx_DAN_MBON(neuron_inx,1,ind);
                Dx_DAN_MBON_Matrix_CSm{neuron_inx,1}(indxx,ind) = Dx_DAN_MBON(neuron_inx,2,ind);
                Dx_DAN_MBON_Matrix_delta{neuron_inx,1}(indxx,ind) = Dx_DAN_MBON(neuron_inx,1,ind) - Dx_DAN_MBON(neuron_inx,2,ind);
            end
        end
        indxx = indxx+1;
    end
    for neuron_inx = 1:6
        MemoryTrace_CSp{neuron_inx,1}(indx,1:2) = Dx_DAN_MBON_Matrix_CSp{neuron_inx,1}(1,1:2);
        MemoryTrace_CSm{neuron_inx,1}(indx,1:2) = Dx_DAN_MBON_Matrix_CSm{neuron_inx,1}(1,1:2);
        MemoryTrace_delta{neuron_inx,1}(indx,1:2) = Dx_DAN_MBON_Matrix_delta{neuron_inx,1}(1,1:2);
        
        MemoryTrace_CSp{neuron_inx,1}(indx,2:indxx) = Dx_DAN_MBON_Matrix_CSp{neuron_inx,1}(:,size(Dx_DAN_MBON,3))';
        MemoryTrace_CSm{neuron_inx,1}(indx,2:indxx) = Dx_DAN_MBON_Matrix_CSm{neuron_inx,1}(:,size(Dx_DAN_MBON,3))';
        MemoryTrace_delta{neuron_inx,1}(indx,2:indxx) = Dx_DAN_MBON_Matrix_delta{neuron_inx,1}(:,size(Dx_DAN_MBON,3))';
    end
    indx = indx+1;
end

%%
paper_size=[6 11];

fig_h=figure('Units','inches', ...
             'Papersize', paper_size, ...
             'PaperPosition', [0 0, paper_size], ...
             'Position',[0.5 0.5 paper_size]);

paper_margin=[(8.5-6.7)/2 0.5];
font_size=7;

% cMap = getPyPlot_cMap('RdBu_r');
J = customcolormap([0 0.25 0.375 0.5 0.625 0.75 1], {'#8B008B','#C71585','#FA8072','#FFFAFA','#48D1CC','#4682B4','#000000'});
% J = customcolormap([0 0.125 0.25 0.375 0.5 0.625 0.75 0.875 1], {'#8B008B','#C71585','#FA8072','#EEE8AA','#FFFAFA','#3CB371','#48D1CC','#4682B4','#000000'});

colormap(J);
axis_h = zeros(6,3);
for neuron_inx = 1:3
    %------------------------------------------------------------------------
    axis_h(1,neuron_inx) = subplot(7,3,neuron_inx);
    imagesc((rest_t_list+135)/3600,odor_valence_list, ...
        MemoryTrace_CSp{neuron_inx,1}(:,2:end),[-40,40])
    
    axis_h(2,neuron_inx) = subplot(7,3,3+neuron_inx);
    imagesc((rest_t_list+135)/3600,odor_valence_list, ...
        MemoryTrace_CSm{neuron_inx,1}(:,2:end),[-40,40])
    
    axis_h(3,neuron_inx) = subplot(7,3,neuron_inx+6);
    imagesc((rest_t_list+135)/3600,odor_valence_list, ...
        MemoryTrace_delta{neuron_inx,1}(:,2:end),[-40,40])
end

for neuron_inx = 4:6
    %------------------------------------------------------------------------
    axis_h(4,neuron_inx-3) = subplot(7,3,neuron_inx+9);
    imagesc((rest_t_list+135)/3600,odor_valence_list, ...
        MemoryTrace_CSp{neuron_inx,1}(:,2:end)-MemoryTrace_CSp{neuron_inx,1}(:,1),[-40 40])
    
    axis_h(5,neuron_inx-3) = subplot(7,3,neuron_inx+12);
    imagesc((rest_t_list+135)/3600,odor_valence_list, ...
        MemoryTrace_CSm{neuron_inx,1}(:,2:end)-MemoryTrace_CSm{neuron_inx,1}(:,1),[-40 40])
        
    axis_h(6,neuron_inx-3) = subplot(7,3,neuron_inx+15);
    imagesc((rest_t_list+135)/3600,odor_valence_list, ...
        MemoryTrace_delta{neuron_inx,1}(:,2:end),[-40 40])
end

set(axis_h,'box','off','YDir','normal');%, ...
    %'YTick',1:2:9,'YTickLabel',-4:2:4, ...
    %'XTick',[1 6.5 12],'XTickLabel',0:0.5:1)

subplot(7,3,1);title('DAN-v1')
subplot(7,3,2);title('DAN-a2')
subplot(7,3,3);title('DAN-a3')

subplot(7,3,13);title('MBON-v1')
subplot(7,3,14);title('MBON-a2')
subplot(7,3,15);title('MBON-a3')

%%
saveas(fig_h, figure_name, 'fig')
print('-dtiff','-r300',[figure_name,'.tif'])
%print('-dpdf',[fullfile(folder_name,out_file_name),'.pdf'])
print('-dmeta',[figure_name,'.emf'])


end
