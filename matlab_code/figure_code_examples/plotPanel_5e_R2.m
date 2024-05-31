clear all
clc
config_path_and_function

is_simplified_model = false;%true;%
%% load parameters
if is_simplified_model
    load('Dx_steady_state_nonlinear_3_27-Mar-2023_2modules') %%Extended Data Fig10j
else
    load('Dx_steady_state_nonlinear_3_27-Mar-2023_3modules') %% Fig5c
end

%%
[exp_para, exp_struct] = experimental_condition([]);

exp_para.session_list = ...
    {'imaging';
    'training';
    'rest';
    'imaging'};


exp_para.session(3).t_ISI = 135; % interval
exp_para.session(3).n = 1; % training cycle % Figure 5d

odor_valence = 0;

for model_inx = 1:2
    switch model_inx
        case 1
            Modeltype = 'Intact Model';
            para_mat_cell = reset_parameter_by_valance(para_mu,mat_lu_cell,odor_valence);
        case 2
            Modeltype = 'Model without Feedback';
            para_mu([12,13,15,18,19],1) = 0; % disconnect MBON-v1 feedbacks 3 module
            para_mat_cell = reset_parameter_by_valance(para_mu,mat_lu_cell,odor_valence);
    end

            
para_mat_cell = reset_parameter_by_valance(para_mu,mat_lu_cell,odor_valence);

training_n_list = 3:3:15;
rest_t_list = (900:900:10800)-135;

[exp_para, exp_struct] = experimental_condition(exp_para);

%% 
Dx_DAN_MBON_Matrix_CSp = cell(6,1);
Dx_DAN_MBON_Matrix_CSm = cell(6,1);
Dx_DAN_MBON_Matrix_delta = cell(6,1);
MemoryTrace_CSp = cell(6,1);
MemoryTrace_CSm = cell(6,1);
MemoryTrace_delta = cell(6,1); 

indx = 1;
for training_n = training_n_list 

    exp_para.session(3).n = training_n;

    indxx = 1;

    for rest_t = rest_t_list

        exp_para.rest_t = rest_t;
        [exp_para, exp_struct] = experimental_condition(exp_para);
        exp_struct(4).t_length = 300;
        [Dx_DAN_MBON,data_all] = Dx_steady_state_MBON_revision(para_mat_cell,exp_struct);

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
figure
sgtitle(Modeltype);

J = customcolormap([0 0.25 0.375 0.5 0.625 0.75 1], {'#8B008B','#C71585','#FA8072','#FFFAFA','#48D1CC','#4682B4','#000000'});

colormap(J);

axis_h = zeros(6,3);
for neuron_inx = 1:3
    %------------------------------------------------------------------------
    axis_h(1,neuron_inx) = subplot(7,3,neuron_inx);
    imagesc((rest_t_list+135)/3600,training_n_list, ...
        MemoryTrace_CSp{neuron_inx,1}(:,2:end),[-40,40])
    
    axis_h(2,neuron_inx) = subplot(7,3,3+neuron_inx);
    imagesc((rest_t_list+135)/3600,training_n_list, ...
        MemoryTrace_CSm{neuron_inx,1}(:,2:end),[-40,40])
    
    axis_h(3,neuron_inx) = subplot(7,3,neuron_inx+6);
    imagesc((rest_t_list+135)/3600,training_n_list, ...
        MemoryTrace_delta{neuron_inx,1}(:,2:end),[-40,40])
end

for neuron_inx = 4:6
    %------------------------------------------------------------------------
    axis_h(4,neuron_inx-3) = subplot(7,3,neuron_inx+9);
    imagesc((rest_t_list+135)/3600,training_n_list, ...
        MemoryTrace_CSp{neuron_inx,1}(:,2:end)-MemoryTrace_CSp{neuron_inx,1}(:,1),[-40 40])
    
    axis_h(5,neuron_inx-3) = subplot(7,3,neuron_inx+12);
    imagesc((rest_t_list+135)/3600,training_n_list, ...
        MemoryTrace_CSm{neuron_inx,1}(:,2:end)-MemoryTrace_CSm{neuron_inx,1}(:,1),[-40 40])
        
    axis_h(6,neuron_inx-3) = subplot(7,3,neuron_inx+15);
    imagesc((rest_t_list+135)/3600,training_n_list, ...
        MemoryTrace_delta{neuron_inx,1}(:,2:end),[-40 40])
end

set(axis_h,'box','off','YDir','normal', ...
    'YTick',-4:2:4,'YTickLabel',-4:2:4, ...
    'XTick',0:1:3,'XTickLabel',0:1:3);

subplot(7,3,1);title('DAN-v1')
subplot(7,3,2);title('DAN-a2')
subplot(7,3,3);title('DAN-a3')

subplot(7,3,13);title('MBON-v1')
subplot(7,3,14);title('MBON-a2')
subplot(7,3,15);title('MBON-a3')
end