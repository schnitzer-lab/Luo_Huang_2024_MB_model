clear all
clc
config_path_and_function
load('Dx_steady_state_nonlinear_3_27-Mar-2023_3modules')

%% Intact model
MemoryTrace_CSp_intact = cell(6,3);
MemoryTrace_CSm_intact = cell(6,3);
MemoryTrace_delta_intact = cell(6,3);

for i = 1:length(para_rand)

    odor_valence = 0;
    para_mat_cell = reset_parameter_by_valance(para_rand(:,i),mat_lu_cell,odor_valence);

    %-----loading model fitted parameters
    [exp_para, exp_struct] = experimental_condition([]);
    [exp_para, exp_struct] = experimental_condition(exp_para);

    exp_para.session_list = ...
        {'imaging';
        'training';
        'rest';
        'imaging'};

    exp_para.session(3).t_ISI = 135; % interval

    training_n_list = 3:3:15;
    rest_t_list = [900-135,10800-135];

    [exp_para, exp_struct] = experimental_condition(exp_para);


    Dx_DAN_MBON_Matrix_CSp = cell(6,1);
    Dx_DAN_MBON_Matrix_CSm = cell(6,1);
    Dx_DAN_MBON_Matrix_delta = cell(6,1);

    indx = 1;
    for training_n = training_n_list

        exp_para.session(3).n = training_n;

        indxx = 1;

        for rest_t = rest_t_list

            exp_para.rest_t = rest_t;
            [exp_para, exp_struct] = experimental_condition(exp_para);
            exp_struct(4).t_length = 300;
            [Dx_DAN_MBON{i},data_all] = Dx_steady_state_MBON_revision(para_mat_cell,exp_struct);

            for neuron_inx = 1:6
                for ind = 1:size(Dx_DAN_MBON{i},3)
                    Dx_DAN_MBON_Matrix_CSp{neuron_inx,1}(indxx,ind) = Dx_DAN_MBON{i}(neuron_inx,1,ind);
                    Dx_DAN_MBON_Matrix_CSm{neuron_inx,1}(indxx,ind) = Dx_DAN_MBON{i}(neuron_inx,2,ind);
                    Dx_DAN_MBON_Matrix_delta{neuron_inx,1}(indxx,ind) = Dx_DAN_MBON{i}(neuron_inx,1,ind) - Dx_DAN_MBON{i}(neuron_inx,2,ind);
                end
            end
            indxx = indxx+1;
        end
        for neuron_inx = 1:6
            MemoryTrace_CSp_intact{neuron_inx,indx}(i,1:2) = Dx_DAN_MBON_Matrix_CSp{neuron_inx,1}(1,1:2);
            MemoryTrace_CSm_intact{neuron_inx,indx}(i,1:2) = Dx_DAN_MBON_Matrix_CSm{neuron_inx,1}(1,1:2);
            MemoryTrace_delta_intact{neuron_inx,indx}(i,1:2) = Dx_DAN_MBON_Matrix_delta{neuron_inx,1}(1,1:2);

            MemoryTrace_CSp_intact{neuron_inx,indx}(i,2:indxx) = Dx_DAN_MBON_Matrix_CSp{neuron_inx,1}(:,size(Dx_DAN_MBON{i},3))';
            MemoryTrace_CSm_intact{neuron_inx,indx}(i,2:indxx) = Dx_DAN_MBON_Matrix_CSm{neuron_inx,1}(:,size(Dx_DAN_MBON{i},3))';
            MemoryTrace_delta_intact{neuron_inx,indx}(i,2:indxx) = Dx_DAN_MBON_Matrix_delta{neuron_inx,1}(:,size(Dx_DAN_MBON{i},3))';
        end
        indx = indx+1;
    end
end

%% No feedback model
MemoryTrace_CSp_noFB = cell(6,3);
MemoryTrace_CSm_noFB = cell(6,3);
MemoryTrace_delta_noFB = cell(6,3);

for i = 1:length(para_rand)

    odor_valence = 0;
    para_rand([12,13,15,18,19],i) = 0; % disconnect MBON-v1 feedbacks
    para_mat_cell = reset_parameter_by_valance(para_rand(:,i),mat_lu_cell,odor_valence);

    %-----loading model fitted parameters
    [exp_para, exp_struct] = experimental_condition([]);
    [exp_para, exp_struct] = experimental_condition(exp_para);

    exp_para.session_list = ...
        {'imaging';
        'training';
        'rest';
        'imaging'};

    exp_para.session(3).t_ISI = 135; % interval

    training_n_list = 3:3:15;
    rest_t_list = [900-135,10800-135];

    [exp_para, exp_struct] = experimental_condition(exp_para);


    Dx_DAN_MBON_Matrix_CSp = cell(6,1);
    Dx_DAN_MBON_Matrix_CSm = cell(6,1);
    Dx_DAN_MBON_Matrix_delta = cell(6,1);

    indx = 1;
    for training_n = training_n_list

        exp_para.session(3).n = training_n;

        indxx = 1;

        for rest_t = rest_t_list

            exp_para.rest_t = rest_t;
            [exp_para, exp_struct] = experimental_condition(exp_para);
            exp_struct(4).t_length = 300;
            [Dx_DAN_MBON{i},data_all] = Dx_steady_state_MBON_revision(para_mat_cell,exp_struct);

            for neuron_inx = 1:6
                for ind = 1:size(Dx_DAN_MBON{i},3)
                    Dx_DAN_MBON_Matrix_CSp{neuron_inx,1}(indxx,ind) = Dx_DAN_MBON{i}(neuron_inx,1,ind);
                    Dx_DAN_MBON_Matrix_CSm{neuron_inx,1}(indxx,ind) = Dx_DAN_MBON{i}(neuron_inx,2,ind);
                    Dx_DAN_MBON_Matrix_delta{neuron_inx,1}(indxx,ind) = Dx_DAN_MBON{i}(neuron_inx,1,ind) - Dx_DAN_MBON{i}(neuron_inx,2,ind);
                end
            end
            indxx = indxx+1;
        end
        for neuron_inx = 1:6
            MemoryTrace_CSp_noFB{neuron_inx,indx}(i,1:2) = Dx_DAN_MBON_Matrix_CSp{neuron_inx,1}(1,1:2);
            MemoryTrace_CSm_noFB{neuron_inx,indx}(i,1:2) = Dx_DAN_MBON_Matrix_CSm{neuron_inx,1}(1,1:2);
            MemoryTrace_delta_noFB{neuron_inx,indx}(i,1:2) = Dx_DAN_MBON_Matrix_delta{neuron_inx,1}(1,1:2);

            MemoryTrace_CSp_noFB{neuron_inx,indx}(i,2:indxx) = Dx_DAN_MBON_Matrix_CSp{neuron_inx,1}(:,size(Dx_DAN_MBON{i},3))';
            MemoryTrace_CSm_noFB{neuron_inx,indx}(i,2:indxx) = Dx_DAN_MBON_Matrix_CSm{neuron_inx,1}(:,size(Dx_DAN_MBON{i},3))';
            MemoryTrace_delta_noFB{neuron_inx,indx}(i,2:indxx) = Dx_DAN_MBON_Matrix_delta{neuron_inx,1}(:,size(Dx_DAN_MBON{i},3))';
        end
        indx = indx+1;
    end
end
%%

for i = 1:5
    for neuron_inx = 1:6
        MemoryTrace_delta_15min{neuron_inx}(1:10000,i)= 1-MemoryTrace_delta_noFB{neuron_inx,i}(:,2)./MemoryTrace_delta_intact{neuron_inx,i}(:,2);
        MemoryTrace_delta_3hr{neuron_inx}(1:10000,i)= 1-MemoryTrace_delta_noFB{neuron_inx,i}(:,3)./MemoryTrace_delta_intact{neuron_inx,i}(:,3);
        A = [];
        A = MemoryTrace_delta_15min{neuron_inx}(1:10000,i);
        A = A(~isnan(A));
        A = A(~isinf(A));
        MemoryTrace_delta_15min_avg(neuron_inx,i) = median(A);
        MemoryTrace_delta_15min_sem(neuron_inx,i) = std(A)/sqrt(12);
        MemoryTrace_delta_15min_std(neuron_inx,i) = std(A);
        
        A = sort(A);
        p = 0.68;
        MemoryTrace_delta_15min_CI{neuron_inx,1}(1,i) = A(round(length(A)*(1+p)/2))-median(A);
        MemoryTrace_delta_15min_CI{neuron_inx,1}(2,i) = median(A) - A(round(length(A)*(1-p)/2));
        
        AA = [];
        AA = MemoryTrace_delta_3hr{neuron_inx}(1:10000,i);
        AA = AA(~isnan(AA));
        AA = AA(~isinf(AA));
        MemoryTrace_delta_3hr_avg(neuron_inx,i) = median(AA);
        MemoryTrace_delta_3hr_sem(neuron_inx,i) = std(AA)/sqrt(12);
        MemoryTrace_delta_3hr_std(neuron_inx,i) = std(AA);
        
        AA = sort(AA);
        p = 0.68;
        MemoryTrace_delta_3hr_CI{neuron_inx,1}(1,i) = AA(round(length(AA)*(1+p)/2))-median(AA);
        MemoryTrace_delta_3hr_CI{neuron_inx,1}(2,i) = median(AA) - AA(round(length(AA)*(1-p)/2));
        
    end
end

%% MBON percentage difference
MBON_v1_diff = [MemoryTrace_delta_15min_avg(4,:)*100, 0, 0, MemoryTrace_delta_3hr_avg(4,:)*100];
MBON_v1_diff(2:3,1:5) = MemoryTrace_delta_15min_CI{4,1}*100;
MBON_v1_diff(2:3,8:12) = MemoryTrace_delta_3hr_CI{4,1}*100;

MBON_a2_diff = [MemoryTrace_delta_15min_avg(5,:)*100, 0, 0, MemoryTrace_delta_3hr_avg(5,:)*100];
MBON_a2_diff(2:3,1:5) = MemoryTrace_delta_15min_CI{5,1}*100;
MBON_a2_diff(2:3,8:12) = MemoryTrace_delta_3hr_CI{5,1}*100;

MBON_a3_diff = [MemoryTrace_delta_15min_avg(6,:)*100, 0, 0, MemoryTrace_delta_3hr_avg(6,:)*100];
MBON_a3_diff(2:3,1:5) = MemoryTrace_delta_15min_CI{6,1}*100;
MBON_a3_diff(2:3,8:12) = MemoryTrace_delta_3hr_CI{6,1}*100;

%% Plot MBONs
for MBON_inx = 1:3
    switch MBON_inx
        case 1
            data = MBON_v1_diff(1,:);
            errhigh = MBON_v1_diff(2,:);
            errlow  = MBON_v1_diff(3,:);
            MBON_name = 'MBON-\gamma1pedc';
        case 2
            data = MBON_a2_diff(1,:);
            errhigh = MBON_a2_diff(2,:);
            errlow  = MBON_a2_diff(3,:);
            MBON_name = 'MBON-\alpha2';
        case 3
            data = MBON_a3_diff(1,:);
            errhigh = MBON_a3_diff(2,:);
            errlow  = MBON_a3_diff(3,:);
            MBON_name = 'MBON-\alpha3';
    end

    subplot(1,3,MBON_inx)
    X = 1:12;
    bar(X,data)
    set(gca,'XTickLabel',[]);
    title(MBON_name,'FontSize', 15)
    box off
    hold on


    er = errorbar(X,data,errlow,errhigh);
    er.Color = [0 0 0];
    er.LineStyle = 'none';
    ylabel('percentage difference with or without feedback','FontSize', 15)

    text(4,-8,'Cycles of training','FontSize', 15)
    text(3,-12,'15-min','FontSize', 15)
    text(10,-12,'3-hr','FontSize', 15)

    text(1,-5,'3','FontSize', 15)
    text(2,-5,'6','FontSize', 15)
    text(3,-5,'9','FontSize', 15)
    text(4,-5,'12','FontSize', 15)
    text(5,-5,'15','FontSize', 15)
    text(8,-5,'3','FontSize', 15)
    text(9,-5,'6','FontSize', 15)
    text(10,-5,'9','FontSize', 15)
    text(11,-5,'12','FontSize', 15)
    text(12,-5,'15','FontSize', 15)

    ylim([0 110])
    hold off

end










