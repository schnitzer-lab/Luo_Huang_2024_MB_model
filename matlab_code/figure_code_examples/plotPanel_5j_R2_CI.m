clear all
clc
config_path_and_function
load('Dx_steady_state_nonlinear_3_27-Mar-2023_3modules')

figure_type_list = {'ACV_EtA','OCT_BEN'};

for figure_type_i= 1:length(figure_type_list)
    figure_type = figure_type_list{figure_type_i};
    switch figure_type
        case 'ACV_EtA'
            odor_valence = [-3.42;-4.58;-5.60]; % ACV/EtA

        case 'OCT_BEN'
            odor_valence = [2.57;2.24;3.35]; % BEN/OCT
    end

    t_extinction_list = [nan 10 120]*60-300-250;

    for i = 1:length(para_rand)

        para_mat_cell = reset_parameter_by_valance(para_rand(:,i),mat_lu_cell,odor_valence);

        shock_strength = 1.5; % shock intensity
        para_mat_cell{3}(1) = para_mat_cell{3}(1)*shock_strength;

        rest_t_list = 3600*3-300-250;

        for t_extinction_i = 1:length(t_extinction_list)
            t_extinction_curr = t_extinction_list(t_extinction_i);
            if isnan(t_extinction_curr)
                % create the defualt experimental condition listexp_para
                [exp_para, exp_struct] = experimental_condition([]);

                exp_para.session_list = ...
                    {'imaging';
                    'training';
                    'imaging';
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
                    'imaging';
                    'rest'
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
                exp_struct(4).t_length = 300;
                exp_struct(4*(exp_para.session(3).n+1)).t_length = 300;

                if exp_para.rest_t >= 0
                    [Dx_DAN_MBON{i},data_all] = Dx_steady_state_MBON_revision(para_mat_cell,exp_struct);
                    % Dx_DAN_MBON: first dimension: DAN,MBON; second dimension: CS+, CS-
                    output_data_matrix(i,:,1:2,:,t_extinction_i) = permute(Dx_DAN_MBON{i}(:,:,:),[4 3 2 5 1]);

                elseif rest_t_list(rest_t_i) <= t_extinction_curr
                    output_data_matrix(i,:,1:2,:,t_extinction_i) = ...
                        output_data_matrix(i,:,1:2,:,1);
                else
                    output_data_matrix(i,:,1:2,:,1) = nan;
                end
            end
            output_data_matrix(i,:,3,:,:) = output_data_matrix(i,:,1,:,:) - output_data_matrix(i,:,2,:,:);
        end
    end

    switch figure_type
        case 'ACV_EtA'
            output_data_matrix_attractive = output_data_matrix; % ACV/EtA

        case 'OCT_BEN'
            output_data_matrix_repulsive = output_data_matrix; % BEN/OCT
    end
end

%%
MemoryTrace_delta_attractive_3hr = cell(1,1);
for i = 1:3
    for neuron_inx = 1:6
        MemoryTrace_CSp_attractive_avg{neuron_inx,1}(i,:) = mean(output_data_matrix_attractive(:,:,1,neuron_inx,i));
        MemoryTrace_CSm_attractive_avg{neuron_inx,1}(i,:) = mean(output_data_matrix_attractive(:,:,2,neuron_inx,i));
        MemoryTrace_delta_attractive_avg{neuron_inx,1}(i,:) = mean(output_data_matrix_attractive(:,:,3,neuron_inx,i));

        MemoryTrace_CSp_repulsive_avg{neuron_inx,1}(i,:) = mean(output_data_matrix_repulsive(:,:,1,neuron_inx,i));
        MemoryTrace_CSm_repulsive_avg{neuron_inx,1}(i,:) = mean(output_data_matrix_repulsive(:,:,2,neuron_inx,i));
        MemoryTrace_delta_repulsive_avg{neuron_inx,1}(i,:) = mean(output_data_matrix_repulsive(:,:,3,neuron_inx,i));

        MemoryTrace_delta_attractive_3hr{neuron_inx,1}(:,i) = output_data_matrix_attractive(:,3,3,neuron_inx,i)./output_data_matrix_attractive(:,2,3,neuron_inx,i);
        MemoryTrace_delta_repulsive_3hr{neuron_inx,1}(:,i) = output_data_matrix_repulsive(:,3,3,neuron_inx,i)./output_data_matrix_repulsive(:,2,3,neuron_inx,i);

        A = MemoryTrace_delta_attractive_3hr{neuron_inx,1}(:,i);
        A = A(~isnan(A));
        A = A(~isinf(A));
        MemoryTrace_delta_attractive_3hr_avg{neuron_inx,1}(:,i) = median(A);

        A = sort(A);
        p = 0.68;
        MemoryTrace_delta_attractive_3hr_CI{neuron_inx,1}(1,i) = A(round(length(A)*(1+p)/2))-median(A);
        MemoryTrace_delta_attractive_3hr_CI{neuron_inx,1}(2,i) = median(A) - A(round(length(A)*(1-p)/2));

        AA = MemoryTrace_delta_repulsive_3hr{neuron_inx,1}(:,i);
        AA = AA(~isnan(AA));
        AA = AA(~isinf(AA));
        MemoryTrace_delta_repulsive_3hr_avg{neuron_inx,1}(:,i) = median(AA);

        AA = sort(AA);
        p = 0.68;
        MemoryTrace_delta_repulsive_3hr_CI{neuron_inx,1}(1,i) = AA(round(length(AA)*(1+p)/2))-median(AA);
        MemoryTrace_delta_repulsive_3hr_CI{neuron_inx,1}(2,i) = median(AA) - AA(round(length(AA)*(1-p)/2));


    end
end

%% plot MBON-a3
MBON_a3_attractive = [MemoryTrace_delta_attractive_3hr_avg{6,1}(1,1), 0, 0, 0, MemoryTrace_delta_attractive_3hr_avg{6,1}(1,2), ...
    0, 0, 0, MemoryTrace_delta_attractive_3hr_avg{6,1}(1,3),0]; 
MBON_a3_repulsive = [0, MemoryTrace_delta_repulsive_3hr_avg{6,1}(1,1), 0, 0, 0, MemoryTrace_delta_repulsive_3hr_avg{6,1}(1,2), ...
    0, 0, 0, MemoryTrace_delta_repulsive_3hr_avg{6,1}(1,3)]; 

MBON_a3_attractive(2:3,:) = [MemoryTrace_delta_attractive_3hr_CI{6,1}(:,1), zeros(2,3), MemoryTrace_delta_attractive_3hr_CI{6,1}(:,2), ...
    zeros(2,3), MemoryTrace_delta_attractive_3hr_CI{6,1}(:,3),zeros(2,1)];
MBON_a3_repulsive(2:3,:) = [zeros(2,1), MemoryTrace_delta_repulsive_3hr_CI{6,1}(:,1), zeros(2,3), MemoryTrace_delta_repulsive_3hr_CI{6,1}(:,2), ...
    zeros(2,3), MemoryTrace_delta_repulsive_3hr_CI{6,1}(:,3)]; 

X = 1:10;
bar(X,MBON_a3_attractive(1,:))
box off
hold on
bar(X,MBON_a3_repulsive(1,:))
er = errorbar(X,MBON_a3_attractive(1,:),MBON_a3_attractive(3,:),MBON_a3_attractive(2,:));
er2 = errorbar(X,MBON_a3_repulsive(1,:),MBON_a3_repulsive(3,:),MBON_a3_repulsive(2,:));

er.Color = [0 0 0];
er.LineStyle = 'none'; 
er2.Color = [0 0 0];
er2.LineStyle = 'none'; 
set(gca,'XTickLabel',[]);
title('Simulated MBON-\alpha3','FontSize', 15)
ylabel('3-hr to 5-min ratio','FontSize', 15)
text(1,1.8,'No Extinction','FontSize', 15)
text(5,1.8,'Extinction @10min','FontSize', 15)
text(9,1.8,'Extinction @2hr','FontSize', 15)

ylim([0 2])
hold off