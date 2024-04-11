function fun_plotPanel_5j(file_parameter,figure_name,figure_type)


%load('/Users/flyspike/Desktop/Fly leanring paper/fly model with CI/version2 - a2 a3 equal tau/figures/Dx_steady_state_nonlinear_3_18-Aug-2021')
%load('figures/Dx_steady_state_nonlinear_3_18-Aug-2021')
load(file_parameter)
addpath('customcolormap\')

%-----loading model fitted parameters
%para_mat_cell = parameter_vec2mat(x_lu,mat_lu_cell);
para_mat_cell_all = parameter_vec2mat(para_mu,mat_lu_cell);
switch figure_type
    case 'ACV_EtA'
        para_mat_cell = para_mat_cell_all;
        para_mat_cell{1} = para_mat_cell{1}(:,1:6);
    case 'OCT_BEN'
        para_mat_cell = para_mat_cell_all;
        para_mat_cell{1} = para_mat_cell{1}(:,[7:9 4:6]);
    otherwise
        error('figure_type should be ACV_EtA and OCT_BEN.')
end

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


%%
t_extinction_list = [nan 10 20 30 60 90 120]*60;
rest_t_list = (300:300:3600*3)-135;
output_data_matrix = nan(length(t_extinction_list),length(rest_t_list),3,6);
% the first two dimensions are row and column of the color maps, 
% the last two dimensions correspond to the row and column of subplots.

for t_extinction_i = 1:length(t_extinction_list)
    t_extinction_curr = t_extinction_list(t_extinction_i);
    if isnan(t_extinction_curr)
        % create the defualt experimental condition listexp_para
        [exp_para, exp_struct] = experimental_condition([]);
        
        exp_para.session_list = ...
            {'imaging';
            'training';
            'rest';
            'imaging'};
        
        exp_para.session(3).t_ISI = 135; % interval
        exp_para.session(3).n = 3; % training cycle 

        exp_para.rest_t = 0;
    else
        [exp_para, exp_struct] = experimental_condition([]);
        
        exp_para.session_list = ...
            {'imaging';
            'training';
            'rest';
            'extinction';
            'rest';
            'imaging'};
        
        exp_para.session(3).t_ISI = 135; % interval
        exp_para.session(3).n = 3; % training cycle 
        
        exp_para.session(4) = exp_para.session(3);
        exp_para.session(4).name = 'extinction';

        exp_para.rest_t = [t_extinction_curr; 0];
        
    end
    


    for rest_t_i = 1:length(rest_t_list)
        if isnan(t_extinction_curr)
            exp_para.rest_t = rest_t_list(rest_t_i);
        else
            exp_para.rest_t(end) = rest_t_list(rest_t_i) - exp_para.rest_t(1) ...
                - exp_para.session(4).n*(exp_para.session(4).t_CS_plus ...
                                        +exp_para.session(4).t_ISI ...
                                        +exp_para.session(4).t_CS_minus ...
                                        +exp_para.session(4).t_ISI);
        end
        
        [exp_para, exp_struct] = experimental_condition(exp_para);
        
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
        if exp_para.rest_t >= 0
            [Dx_DAN_MBON,data_all] = Dx_steady_state_MBON_0301_2023(para_mat_cell,exp_struct);
            % Dx_DAN_MBON: first dimension: DAN,MBON; second dimension: CS+, CS-
            output_data_matrix(t_extinction_i,rest_t_i,1:2,:) = permute(Dx_DAN_MBON(:,:,end),[4 3 2 1]);
        elseif rest_t_list(rest_t_i) <= t_extinction_curr
            output_data_matrix(t_extinction_i,rest_t_i,1:2,:) = ...
                output_data_matrix(1,rest_t_i,1:2,:);
        else
            output_data_matrix(t_extinction_i,rest_t_i,1:2,:) = nan;
        end
    end
    output_data_matrix(:,:,3,:) = output_data_matrix(:,:,1,:) - output_data_matrix(:,:,2,:);
    
end

%%
paper_size=[11 6];

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
axis_h = zeros(3,6);
for row_i = 1:3
    for col_i = 1:6
        axis_h(row_i,col_i) = subplot(3,6,(row_i-1)*6+col_i); hold on;
        imagesc(rest_t_list/3600,1:length(t_extinction_list), ...
                output_data_matrix(:,:,row_i,col_i),[-40,40])
    
    end
end


set(axis_h,'box','off','YDir','normal', ...
    'YLim',[0.5,length(t_extinction_list) + 0.5], ...
    'YTick',1:length(t_extinction_list), ...
    'YTickLabel',[{'Ctr'};cellstr(num2str(t_extinction_list(2:end)'/60))], ...
    'XTick',0:3)

title(axis_h(1,1),'DAN-v1')
title(axis_h(1,2),'DAN-a2')
title(axis_h(1,3),'DAN-a3')

title(axis_h(1,4),'MBON-v1')
title(axis_h(1,5),'MBON-a2')
title(axis_h(1,6),'MBON-a3')

%%

if ~isempty(figure_name)
    saveas(fig_h, figure_name, 'fig')
    print('-dtiff','-r300',[figure_name,'.tif'])
    %print('-dpdf',[fullfile(folder_name,out_file_name),'.pdf'])
    print('-dmeta',[figure_name,'.emf'])
end

end
